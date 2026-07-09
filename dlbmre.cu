// dlb_compact_mre.cu
// nvcc -O3 -std=c++17 -arch=sm_70 dlb_compact_mre.cu -o dlb_test
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/copy.h>
#include <thrust/scan.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>
#include <chrono>
#include <cassert>

#ifndef DLB_BLOCK_THREADS
#define DLB_BLOCK_THREADS 256
#endif

// ---------- 小工具 ----------
#define CUDA_CHECK(x) do{cudaError_t err=(x); if(err!=cudaSuccess){ \
  fprintf(stderr,"CUDA error %s @ %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
  std::abort(); }}while(0)

inline int divCeil(int a, int b){ return (a + b - 1) / b; }

struct GPUTimer {
  cudaEvent_t a{}, b{};
  GPUTimer(){ cudaEventCreate(&a); cudaEventCreate(&b); }
  ~GPUTimer(){ cudaEventDestroy(a); cudaEventDestroy(b); }
  void start(cudaStream_t s=0){ cudaEventRecord(a, s); }
  float stop(cudaStream_t s=0){ cudaEventRecord(b, s); cudaEventSynchronize(b);
    float ms=0; cudaEventElapsedTime(&ms,a,b); return ms; }
};

inline uint64_t fnv1a64(const void* data, size_t nbytes){
  const uint8_t* p = static_cast<const uint8_t*>(data);
  uint64_t h = 1469598103934665603ull;
  for(size_t i=0;i<nbytes;i++){ h ^= p[i]; h *= 1099511628211ull; }
  return h;
}

// 有些环境没有 __nanosleep，给个降级版本
__device__ __forceinline__ void nanosleep_fallback(unsigned int cycles){
#if defined(__CUDA_ARCH__)
  #if (__CUDACC_VER_MAJOR__ >= 11)
    __nanosleep(cycles);
  #else
    // busy wait
    for (unsigned int i=0;i<cycles;i++){ __asm__ __volatile__("nop;"); }
  #endif
#else
  (void)cycles;
#endif
}

// ---------- 块内 exclusive 扫描 ----------
template<int NT>
__device__ __forceinline__ int blockExclusiveScan(int v, volatile int* sh){
  int tid = threadIdx.x;
  sh[tid] = v;
  __syncthreads();
  // upsweep (inclusive)
  for(int off=1; off<NT; off<<=1){
    int t = (tid>=off) ? sh[tid-off] : 0;
    __syncthreads();
    sh[tid] += t;
    __syncthreads();
  }
  int total = sh[NT-1];
  int excl  = (tid==0)?0:sh[tid-1];
  __syncthreads();
  sh[tid] = excl;
  __syncthreads();
  return total;
}

// ---------- DLB 单核扫描+压缩 ----------
extern "C" __global__ void compactDecoupledLookback(
    const int N,
    const int* __restrict__ in,   // 含无效项：in[i] >= 0 表示保留
    int* __restrict__ out,        // 紧凑输出
    int* __restrict__ blockSums,  // [numBlocks]
    int* __restrict__ blockStatus,// [numBlocks] 0/1
    int* __restrict__ aliveCount  // [1]
){
  const int gtid = blockIdx.x * blockDim.x + threadIdx.x;
  __shared__ int sh[DLB_BLOCK_THREADS];

  int flag=0, val=-1;
  if (gtid < N){ val = in[gtid]; flag = (val >= 0); }

  int bTotal = blockExclusiveScan<DLB_BLOCK_THREADS>(flag, sh);
  int tExcl  = sh[threadIdx.x];

  if (threadIdx.x==0){
    blockSums[blockIdx.x] = bTotal;
    __threadfence();
    blockStatus[blockIdx.x] = 1;
    __threadfence();
  }
  __syncthreads();

  int bOffset = 0;
  int look = 1;
  int idx = blockIdx.x - look;
  while (idx >= 0){
    while (atomicAdd(&blockStatus[idx], 0) == 0) { nanosleep_fallback(64); }
    bOffset += blockSums[idx];
    look <<= 1;
    idx = blockIdx.x - look;
  }
  __syncthreads();

  if (gtid < N && flag){
    int outIdx = bOffset + tExcl;
    out[outIdx] = val;
  }

  if (gtid == N - 1){
    int alive = bOffset + bTotal;
    *aliveCount = alive;
  }
}

// ---------- 基准：生成随机 nodeIndices（部分为负，代表被移除） ----------
void make_random_input(std::vector<int>& h, int N, float keep_ratio, uint32_t seed=42){
  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> U(0.0f, 1.0f);
  h.resize(N);
  for (int i=0;i<N;i++){
    if (U(rng) < keep_ratio) h[i] = i; // 有效：用 i 作为值
    else h[i] = -1;                     // 无效
  }
}

// ---------- Baseline: thrust 扫描 + copy_if ----------
float run_baseline(const thrust::device_vector<int>& dIn,
                   thrust::device_vector<int>& dOut,
                   int& alive)
{
  GPUTimer t;
  const int N = (int)dIn.size();

  // flags
  thrust::device_vector<int> dFlag(N);
  t.start();
  thrust::transform(dIn.begin(), dIn.end(), dFlag.begin(),
    [] __device__ (int x){ return x>=0 ? 1 : 0; });
  float t_flag = t.stop();

  // exclusive scan
  thrust::device_vector<int> dScan(N);
  t.start();
  thrust::exclusive_scan(dFlag.begin(), dFlag.end(), dScan.begin());
  float t_scan = t.stop();

  // total alive
  int lastFlag=0, lastScan=0;
  CUDA_CHECK(cudaMemcpy(&lastFlag, thrust::raw_pointer_cast(dFlag.data()+N-1), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(&lastScan, thrust::raw_pointer_cast(dScan.data()+N-1), sizeof(int), cudaMemcpyDeviceToHost));
  alive = lastScan + lastFlag;

  // scatter (copy_if 也可)
  t.start();
  // 用 zip 简单搬运（保持稳定）
  thrust::for_each_n(thrust::device, thrust::make_counting_iterator<int>(0), N,
    [in=thrust::raw_pointer_cast(dIn.data()),
     scan=thrust::raw_pointer_cast(dScan.data()),
     out=thrust::raw_pointer_cast(dOut.data()),
     N] __device__ (int i){
      int v = in[i];
      if (v>=0){
        int pos = scan[i];
        out[pos] = v;
      }
    });
  float t_scatter = t.stop();

  return t_flag + t_scan + t_scatter; // 作为 baseline 的“压缩阶段时间”
}

// ---------- DLB: 单核 ----------
float run_dlb(const thrust::device_vector<int>& dIn,
              thrust::device_vector<int>& dOut,
              int& alive)
{
  GPUTimer t;
  const int N = (int)dIn.size();
  int nBlocks = divCeil(N, DLB_BLOCK_THREADS);

  thrust::device_vector<int> dBlockSums(nBlocks, 0);
  thrust::device_vector<int> dBlockStatus(nBlocks, 0);
  thrust::device_vector<int> dAlive(1, 0);

  void* args[] = {
    (void*)&N,
    (void*)thrust::raw_pointer_cast(dIn.data()),
    (void*)thrust::raw_pointer_cast(dOut.data()),
    (void*)thrust::raw_pointer_cast(dBlockSums.data()),
    (void*)thrust::raw_pointer_cast(dBlockStatus.data()),
    (void*)thrust::raw_pointer_cast(dAlive.data())
  };

  t.start();
  CUDA_CHECK(cudaLaunchKernel((void*)compactDecoupledLookback,
    dim3(nBlocks), dim3(DLB_BLOCK_THREADS), nullptr, 0, args));
  CUDA_CHECK(cudaDeviceSynchronize());
  float ms = t.stop();

  CUDA_CHECK(cudaMemcpy(&alive, thrust::raw_pointer_cast(dAlive.data()),
    sizeof(int), cudaMemcpyDeviceToHost));
  return ms;
}

// ---------- 校验/哈希 ----------
uint64_t hash_vec(const thrust::device_vector<int>& d, int n){
  thrust::host_vector<int> h(n);
  CUDA_CHECK(cudaMemcpy(h.data(), thrust::raw_pointer_cast(d.data()),
    n*sizeof(int), cudaMemcpyDeviceToHost));
  return fnv1a64(h.data(), n*sizeof(int));
}

int main(int argc, char** argv){
  int N = 1<<22;          // 元素数（模拟“代表”数）
  float keep = 0.5f;      // 存活比例（合并后剩余的代表比例）
  if (argc>=2) N = std::atoi(argv[1]);
  if (argc>=3) keep = std::atof(argv[2]);

  printf("N = %d, keep_ratio = %.2f, block = %d\n", N, keep, DLB_BLOCK_THREADS);

  // 生成输入
  std::vector<int> hIn;
  make_random_input(hIn, N, keep, /*seed*/123);

  thrust::device_vector<int> dIn = hIn;
  thrust::device_vector<int> dOutBase(N, -2);
  thrust::device_vector<int> dOutDLB (N, -3);

  // Baseline
  int aliveBase=0;
  float msBase = run_baseline(dIn, dOutBase, aliveBase);
  uint64_t hBase = hash_vec(dOutBase, aliveBase);

  // DLB
  int aliveDLB=0;
  float msDLB = run_dlb(dIn, dOutDLB, aliveDLB);
  uint64_t hDLB = hash_vec(dOutDLB, aliveDLB);

  // 打印
  printf("Baseline (scan+scatter): %.3f ms, alive = %d, hash = %016llx\n",
    msBase, aliveBase, (unsigned long long)hBase);
  printf("DLB (1 kernel)         : %.3f ms, alive = %d, hash = %016llx\n",
    msDLB, aliveDLB, (unsigned long long)hDLB);

  bool ok = (aliveBase==aliveDLB) && (hBase==hDLB);
  printf("CONSISTENCY: %s\n", ok? "OK ✅" : "MISMATCH ❌");

  return ok? 0 : 1;
}
