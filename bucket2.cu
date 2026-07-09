// morton_sort_benchmark_debug.cu
// nvcc -O2 -std=c++17 morton_sort_benchmark_debug.cu -o buckettest_dbg
//
// 运行示例：
//   ./buckettest_dbg                        // 默认 N=100000, BINS=128, streams=1, dumpN=0, verbose=1
//   ./buckettest_dbg 1000000 1024 4 16 1    // 指定 N, BINS, streams, dump前16项, verbose=1
//
// 说明：
//   方法A：全局一次 radix sort (Thrust/CUB)
//   方法B：粗分桶(取高位) + 桶内 radix sort
//   本调试版会在每一步输出“阶段开始/结束/耗时/显存”等，帮助定位卡住位置。

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/copy.h>



#define CUDA_CHECK(expr) do { \
  cudaError_t __err = (expr); \
  if (__err != cudaSuccess) { \
    fprintf(stderr, "\n[CUDA ERROR] %s failed at %s:%d : %s\n", \
      #expr, __FILE__, __LINE__, cudaGetErrorString(__err)); \
    fflush(stderr); std::exit(1); \
  } \
} while(0)

static void print_mem(const char* tag) {
  size_t free_b=0,total_b=0;
  cudaMemGetInfo(&free_b,&total_b);
  printf("[MEM] %s: free=%.2f MB / total=%.2f MB\n",
         tag, free_b/1048576.0, total_b/1048576.0);
  fflush(stdout);
}

static void print_dev_info() {
  int dev=0; CUDA_CHECK(cudaGetDevice(&dev));
  cudaDeviceProp p{}; CUDA_CHECK(cudaGetDeviceProperties(&p,dev));
  printf("[GPU] %s, CC %d.%d, SMs=%d, GlobalMem=%.1f GB\n",
         p.name, p.major, p.minor, p.multiProcessorCount, p.totalGlobalMem/1073741824.0);
  fflush(stdout);
}

// ------- Morton(30-bit, 每轴10bit) -------
__host__ __device__ inline uint32_t part1by2(uint32_t x) {
  x &= 0x000003ffu;
  x = (x ^ (x << 16)) & 0x030000FFu;
  x = (x ^ (x << 8 )) & 0x0300F00Fu;
  x = (x ^ (x << 4 )) & 0x030C30C3u;
  x = (x ^ (x << 2 )) & 0x09249249u;
  return x;
}
__host__ __device__ inline uint32_t morton3D_10bits(uint32_t xi,uint32_t yi,uint32_t zi){
  return (part1by2(zi)<<2) | (part1by2(yi)<<1) | (part1by2(xi));
}
__host__ __device__ inline uint32_t quant10(float v,float vmin,float vmax){
  float t=(v - vmin)/(vmax - vmin);
  t = fminf(fmaxf(t,0.0f),1.0f);
  return (uint32_t)(t*1023.0f + 0.5f);
}

// ------- Kernels -------
__global__ void compute_bin_ids(const uint32_t* __restrict__ keys,
                                uint32_t* __restrict__ bin_ids,
                                uint32_t N, int topBits) {
  uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i>=N) return;
  uint32_t k = keys[i];
  bin_ids[i] = (topBits==0)? 0u : (k >> (32 - topBits));
}
__global__ void histogram_bins(const uint32_t* __restrict__ bin_ids,
                               uint32_t* __restrict__ hist, uint32_t N) {
  uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i>=N) return;
  atomicAdd(hist + bin_ids[i], 1u);
}
__global__ void scatter_by_bins(const uint32_t* __restrict__ bin_ids,
                                const uint32_t* __restrict__ in_keys,
                                const uint32_t* __restrict__ in_vals,
                                uint32_t* __restrict__ out_keys,
                                uint32_t* __restrict__ out_vals,
                                uint32_t* __restrict__ bin_write_ptrs,
                                uint32_t N) {
  uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i>=N) return;
  uint32_t b = bin_ids[i];
  uint32_t pos = atomicAdd(bin_write_ptrs + b, 1u);
  out_keys[pos] = in_keys[i];
  out_vals[pos] = in_vals[i];
}

static inline bool is_power_of_two(uint32_t x){ return x && ((x&(x-1))==0); }
static int ilog2_uint(uint32_t x){ int r=-1; while(x){ x>>=1; ++r; } return r; }

// 生成随机 Morton 数据
static void gen_data(thrust::host_vector<uint32_t>& hkeys,
                     thrust::host_vector<uint32_t>& hvals,
                     uint32_t N, uint64_t seed=20250917ULL) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<float> uni(0.0f,1.0f);
  const float xmin=-100.f,xmax=300.f, ymin=-100.f,ymax=300.f, zmin=-100.f,zmax=300.f;
  hkeys.resize(N); hvals.resize(N);
  for(uint32_t i=0;i<N;++i){
    float x=xmin+(xmax-xmin)*uni(rng);  //归一化
    float y=ymin+(ymax-ymin)*uni(rng);
    float z=zmin+(zmax-zmin)*uni(rng);
    uint32_t xi=quant10(x,xmin,xmax);
    uint32_t yi=quant10(y,ymin,ymax);
    uint32_t zi=quant10(z,zmin,zmax);
    hkeys[i]=morton3D_10bits(xi,yi,zi);
    hvals[i]=i;
  }
}

static bool check_equal(const thrust::host_vector<uint32_t>& ak,
                        const thrust::host_vector<uint32_t>& av,
                        const thrust::host_vector<uint32_t>& bk,
                        const thrust::host_vector<uint32_t>& bv,
                        size_t &mismatch_idx) {
  if (ak.size()!=bk.size() || av.size()!=bv.size()) { mismatch_idx=(size_t)-1; return false; }
  size_t N=ak.size();
  for(size_t i=0;i<N;++i){
    if (ak[i]!=bk[i] || av[i]!=bv[i]) { mismatch_idx=i; return false; }
  }
  return true;
}

int main(int argc,char**argv){
  uint32_t N       = (argc>=2)? (uint32_t)std::stoul(argv[1]) : 100000u;
  uint32_t BINS    = (argc>=3)? (uint32_t)std::stoul(argv[2]) : 128u;
  int      streams = (argc>=4)? std::max(1, std::stoi(argv[3])) : 1;
  int      dumpN   = (argc>=5)? std::max(0, std::stoi(argv[4])) : 0; // dump 前N个
  int      verbose = (argc>=6)? std::stoi(argv[5]) : 1;             // 1=详细日志
 int devCount = 0;
cudaGetDeviceCount(&devCount);
if (devCount == 0) {
    printf("未检测到可用的 NVIDIA GPU！\n");
    return 1;
}
// 直接指定第0个 NVIDIA GPU
cudaSetDevice(0);

  if (!is_power_of_two(BINS) || BINS>(1u<<16)){
    printf("BINS 必须为2的幂且<=65536，当前=%u\n",BINS); fflush(stdout); return 1;
  }
  int topBits = ilog2_uint(BINS);

  printf("=== buckettest_dbg ===\n");
  print_dev_info();
  print_mem("initial");
  printf("Config: N=%u, BINS=%u (topBits=%d), streams=%d, dumpN=%d, verbose=%d\n",
         N,BINS,topBits,streams,dumpN,verbose);
  fflush(stdout);

  cudaEvent_t t0,t1; CUDA_CHECK(cudaEventCreate(&t0)); CUDA_CHECK(cudaEventCreate(&t1));
  float ms=0.0f;

  // --- 阶段0：生成数据 ---
  printf("[Stage 0] 生成随机 Morton 数据...\n"); fflush(stdout);
  thrust::host_vector<uint32_t> hkeys, hvals;
  CUDA_CHECK(cudaEventRecord(t0));
  gen_data(hkeys,hvals,N);
  CUDA_CHECK(cudaEventRecord(t1));
  CUDA_CHECK(cudaEventSynchronize(t1));
  CUDA_CHECK(cudaEventElapsedTime(&ms,t0,t1));
  printf("[Stage 0] 完成，用时=%.3f ms\n", ms); fflush(stdout);

  if (verbose && dumpN>0){
    printf("  示例 key[0:%d):", dumpN);
    for(int i=0;i<dumpN && i<(int)N;++i) printf(" %08x", hkeys[i]);
    printf("\n"); fflush(stdout);
  }

  // --- 阶段1：拷贝到GPU ---
  printf("[Stage 1] 拷贝数据到 GPU...\n"); fflush(stdout);
  CUDA_CHECK(cudaEventRecord(t0));
  thrust::device_vector<uint32_t> d_keys=hkeys, d_vals=hvals;
  CUDA_CHECK(cudaEventRecord(t1));
  CUDA_CHECK(cudaEventSynchronize(t1));
  CUDA_CHECK(cudaEventElapsedTime(&ms,t0,t1));
  print_mem("after H2D");
  printf("[Stage 1] 完成，用时=%.3f ms\n", ms); fflush(stdout);

  // --- 方法A：全局一次sort ---
  printf("[Stage 2A] 全局 radix sort (A) 开始...\n"); fflush(stdout);
  thrust::device_vector<uint32_t> A_keys=d_keys, A_vals=d_vals;
  CUDA_CHECK(cudaEventRecord(t0));
  thrust::sort_by_key(A_keys.begin(), A_keys.end(), A_vals.begin());
  CUDA_CHECK(cudaEventRecord(t1));
  CUDA_CHECK(cudaEventSynchronize(t1));
  CUDA_CHECK(cudaEventElapsedTime(&ms,t0,t1));
  float tA = ms;
  printf("[Stage 2A] 完成，用时=%.3f ms\n", tA); fflush(stdout);

  // --- 方法B：分桶+桶内排序 ---
  printf("[Stage 2B] 分桶 + 桶内排序 (B) 准备内存...\n"); fflush(stdout);
  thrust::device_vector<uint32_t> B_in_keys=d_keys, B_in_vals=d_vals;
  thrust::device_vector<uint32_t> B_bin_ids(N);
  thrust::device_vector<uint32_t> B_hist(BINS,0u);
  thrust::device_vector<uint32_t> B_offsets(BINS,0u);
  thrust::device_vector<uint32_t> B_write_ptrs(BINS,0u);
  thrust::device_vector<uint32_t> B_grouped_keys(N), B_grouped_vals(N);
  print_mem("before Stage 2B");

  const int block=256;
  const int gridN=(N+block-1)/block;

  cudaEvent_t tB0,tBmid,tB1;
  CUDA_CHECK(cudaEventCreate(&tB0));
  CUDA_CHECK(cudaEventCreate(&tBmid));
  CUDA_CHECK(cudaEventCreate(&tB1));

  printf("[Stage 2B] 开始分桶阶段 (计算桶号/直方图/前缀和/Scatter)...\n"); fflush(stdout);
  CUDA_CHECK(cudaEventRecord(tB0));

  // 2B-1 计算桶号
  printf("  [2B-1] compute_bin_ids<<<%d,%d>>> ...\n", gridN, block); fflush(stdout);
  compute_bin_ids<<<gridN,block>>>(thrust::raw_pointer_cast(B_in_keys.data()),
                                   thrust::raw_pointer_cast(B_bin_ids.data()),
                                   N, topBits);
  CUDA_CHECK(cudaGetLastError()); CUDA_CHECK(cudaDeviceSynchronize());

  // 2B-2 直方图
  printf("  [2B-2] histogram_bins<<<%d,%d>>> ...\n", gridN, block); fflush(stdout);
  histogram_bins<<<gridN,block>>>(thrust::raw_pointer_cast(B_bin_ids.data()),
                                  thrust::raw_pointer_cast(B_hist.data()), N);
  CUDA_CHECK(cudaGetLastError()); CUDA_CHECK(cudaDeviceSynchronize());

  // 2B-3 前缀和 -> offsets
  printf("  [2B-3] exclusive_scan (B_hist -> B_offsets) ...\n"); fflush(stdout);
  thrust::exclusive_scan(B_hist.begin(), B_hist.end(), B_offsets.begin());
  CUDA_CHECK(cudaDeviceSynchronize());

  // 2B-4 init write_ptrs = offsets
  printf("  [2B-4] copy offsets -> write_ptrs ...\n"); fflush(stdout);
  thrust::copy(B_offsets.begin(), B_offsets.end(), B_write_ptrs.begin());
  CUDA_CHECK(cudaDeviceSynchronize());

  // 2B-5 Scatter
  printf("  [2B-5] scatter_by_bins<<<%d,%d>>> ...\n", gridN, block); fflush(stdout);
  scatter_by_bins<<<gridN,block>>>(
      thrust::raw_pointer_cast(B_bin_ids.data()),
      thrust::raw_pointer_cast(B_in_keys.data()),
      thrust::raw_pointer_cast(B_in_vals.data()),
      thrust::raw_pointer_cast(B_grouped_keys.data()),
      thrust::raw_pointer_cast(B_grouped_vals.data()),
      thrust::raw_pointer_cast(B_write_ptrs.data()),
      N);
  CUDA_CHECK(cudaGetLastError()); CUDA_CHECK(cudaDeviceSynchronize());

  CUDA_CHECK(cudaEventRecord(tBmid));
  CUDA_CHECK(cudaEventSynchronize(tBmid));
  float tB_bucket_ms=0.0f;
  CUDA_CHECK(cudaEventElapsedTime(&tB_bucket_ms,tB0,tBmid));
  printf("[Stage 2B] 分桶阶段完成，用时=%.3f ms\n", tB_bucket_ms); fflush(stdout);

  // 桶信息拉回 host，方便日志
  thrust::host_vector<uint32_t> h_offsets = B_offsets;
  thrust::host_vector<uint32_t> h_counts(BINS);
  for (uint32_t i=0;i<BINS;++i){
    uint32_t start = h_offsets[i];
    uint32_t end = (i+1<BINS)? h_offsets[i+1] : N;
    h_counts[i] = end - start;
  }
  if (verbose){
    uint64_t maxc=0, minc=UINT64_MAX, nonempty=0;
    for(uint32_t i=0;i<BINS;++i){
      if (h_counts[i]) { nonempty++; maxc=std::max<uint64_t>(maxc,h_counts[i]); minc=std::min<uint64_t>(minc,h_counts[i]); }
    }
    double avg = (double)N / (double)BINS;
    printf("  桶统计: 非空桶=%u/%u, 平均=%.1f, 最小=%llu, 最大=%llu\n",
           (unsigned)nonempty,(unsigned)BINS, avg,
           (unsigned long long) (nonempty?minc:0),
           (unsigned long long) maxc);
    fflush(stdout);
  }

  // 2B-6 桶内排序
  printf("[Stage 2B] 桶内排序开始 (streams=%d)...\n", streams); fflush(stdout);
  std::vector<cudaStream_t> ss(streams);
  for(int s=0;s<streams;++s) CUDA_CHECK(cudaStreamCreate(&ss[s]));

  CUDA_CHECK(cudaEventRecord(t0));
  for (uint32_t b=0;b<BINS;++b){
    uint32_t start = h_offsets[b];
    uint32_t cnt   = h_counts[b];
    if (cnt<=1) continue;
    int sid = (streams>0)? (b%streams) : 0;
    auto kb = B_grouped_keys.begin() + start;
    auto ke = kb + cnt;
    auto vb = B_grouped_vals.begin() + start;
    // 分桶内调用 thrust 排序（绑定流）
    thrust::sort_by_key(thrust::cuda::par.on(ss[sid]), kb, ke, vb);
  }
  for (int s=0;s<streams;++s) CUDA_CHECK(cudaStreamSynchronize(ss[s]));
  CUDA_CHECK(cudaEventRecord(t1));
  CUDA_CHECK(cudaEventSynchronize(t1));
  float tB_sort_ms=0.0f;
  CUDA_CHECK(cudaEventElapsedTime(&tB_sort_ms,t0,t1));

  float tB_total_ms=0.0f;
  CUDA_CHECK(cudaEventRecord(tB1));
  CUDA_CHECK(cudaEventSynchronize(tB1));
  CUDA_CHECK(cudaEventElapsedTime(&tB_total_ms,tB0,tB1));

  for (int s=0;s<streams;++s) CUDA_CHECK(cudaStreamDestroy(ss[s]));

  printf("[Stage 2B] 桶内排序完成，用时=%.3f ms; 总计(B)=%.3f ms\n",
         tB_sort_ms, tB_total_ms); fflush(stdout);

  // --- 校验 ---
  printf("[Stage 3] 校验与打印结果...\n"); fflush(stdout);
  thrust::host_vector<uint32_t> hA_keys=A_keys, hA_vals=A_vals;
  thrust::host_vector<uint32_t> hB_keys=B_grouped_keys, hB_vals=B_grouped_vals;

  size_t mismatch=0;
  bool ok = check_equal(hA_keys,hA_vals,hB_keys,hB_vals,mismatch);
  printf("\n===== Benchmark Result (DEBUG) =====\n");
  printf("Method A (Global radix sort): %.3f ms\n", tA);
  printf("Method B (Coarse bin + per-bin sort): total=%.3f ms (binning=%.3f ms, per-bin-sort=%.3f ms)\n",
         tB_total_ms, tB_bucket_ms, tB_sort_ms);
  if (tB_total_ms>0.0f) printf("Speedup (A/B): %.3f x\n", tA / tB_total_ms);
  printf("Correctness: %s\n", ok? "PASS":"FAIL");
  if (!ok){
    printf("First mismatch at idx=%zu\n", mismatch);
    if (verbose){
      size_t i=mismatch;
      printf("  A: key=%08x val=%u\n", hA_keys[i], hA_vals[i]);
      printf("  B: key=%08x val=%u\n", hB_keys[i], hB_vals[i]);
    }
  }
  if (verbose && dumpN>0){
    int K = std::min<int>(dumpN, (int)N);
    printf("\n  A_keys[0:%d):", K);
    for(int i=0;i<K;++i) printf(" %08x", hA_keys[i]); printf("\n");
    printf("  B_keys[0:%d):", K);
    for(int i=0;i<K;++i) printf(" %08x", hB_keys[i]); printf("\n");
  }
  fflush(stdout);

  return ok? 0:2;
}
