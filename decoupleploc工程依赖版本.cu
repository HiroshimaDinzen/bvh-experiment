//段0 ======= utils_timing_hash.hpp =======
#pragma once
#include <cuda_runtime.h>
#include <cstdint>
#include <vector>
#include <cstdio>
#include <cstring>

#define CUDA_CHECK(x) do { cudaError_t err=(x); if (err!=cudaSuccess){ \
  fprintf(stderr,"CUDA error %s @ %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
  std::abort(); } } while(0)

struct CudaTimer {
  cudaEvent_t beg{}, end{};
  CudaTimer()  { cudaEventCreate(&beg); cudaEventCreate(&end); }
  ~CudaTimer() { cudaEventDestroy(beg); cudaEventDestroy(end); }
  void start(cudaStream_t s=0){ cudaEventRecord(beg, s); }
  float stop(cudaStream_t s=0){ cudaEventRecord(end, s); cudaEventSynchronize(end);
    float ms=0; cudaEventElapsedTime(&ms,beg,end); return ms; }
};

// 64-bit FNV-1a 哈希，跨平台稳定（用于一致性校验）
inline uint64_t fnv1a64(const void* data, size_t nbytes){
  const uint8_t* p = static_cast<const uint8_t*>(data);
  uint64_t h = 1469598103934665603ull;
  for(size_t i=0;i<nbytes;i++){ h ^= p[i]; h *= 1099511628211ull; }
  return h;
}

// 从 device 复制到 host（int/float 等简单类型）
template<typename T>
inline std::vector<T> d2h_copy(const T* dptr, size_t n){
  std::vector<T> h(n);
  CUDA_CHECK(cudaMemcpy(h.data(), dptr, n*sizeof(T), cudaMemcpyDeviceToHost));
  return h;
}

//段1
// ======= compact_dlb.cuh =======
#pragma once
#include <cuda_runtime.h>

#ifndef DLB_BLOCK_THREADS
#define DLB_BLOCK_THREADS 256
#endif

template<int NT>
__device__ __forceinline__ int blockExclusiveScan(int v, volatile int* shmem){
  int tid = threadIdx.x;
  shmem[tid] = v;
  __syncthreads();
  for(int off=1; off<NT; off<<=1){
    int t = (tid>=off) ? shmem[tid-off] : 0;
    __syncthreads();
    shmem[tid] += t;
    __syncthreads();
  }
  int total = shmem[NT-1];
  int excl  = (tid==0)?0:shmem[tid-1];
  __syncthreads();
  shmem[tid] = excl;
  __syncthreads();
  return total;
}

extern "C" __global__ void compactDecoupledLookback(
  const int numberOfClusters,
  const int* __restrict__ nodeIndicesIn,
  int* __restrict__ nodeIndicesOut,
  int* __restrict__ blockSums,
  int* __restrict__ blockStatus,
  int* __restrict__ newAliveCount
){
  const int gtid = blockIdx.x * blockDim.x + threadIdx.x;
  __shared__ int shmem[DLB_BLOCK_THREADS];

  int flag=0, val=-1;
  if (gtid < numberOfClusters){
    val  = nodeIndicesIn[gtid];
    flag = (val >= 0);
  }

  int blockTotal = blockExclusiveScan<DLB_BLOCK_THREADS>(flag, shmem);
  int threadExcl = shmem[threadIdx.x];

  if (threadIdx.x==0){
    blockSums[blockIdx.x] = blockTotal;
    __threadfence();
    blockStatus[blockIdx.x] = 1;
    __threadfence();
  }
  __syncthreads();

  int blockOffset = 0;
  int lookback = 1;
  int idx = blockIdx.x - lookback;
  while (idx >= 0){
    while (atomicAdd(&blockStatus[idx], 0) == 0){
      __nanosleep(64);
    }
    blockOffset += blockSums[idx];
    lookback <<= 1;
    idx = blockIdx.x - lookback;
  }
  __syncthreads();

  if (gtid < numberOfClusters && flag){
    int outIdx = blockOffset + threadExcl;
    nodeIndicesOut[outIdx] = val;
  }

  if (gtid == numberOfClusters - 1){
    int alive = blockOffset + blockTotal;
    *newAliveCount = alive;
  }
}
// 段2======= bench_structs.hpp =======
#pragma once
#include <cstdint>
#include <vector>
#include <string>

struct KernelTimes {
  float genNbh = 0.f;
  float merge  = 0.f;
  float localScan = 0.f;
  float globalScan = 0.f;
  float compact = 0.f;   // baseline: 三核总和；DLB: 单核时间
  float total() const { return genNbh + merge + localScan + globalScan + compact; }
};

struct HashSummary {
  uint64_t hParents = 0, hLeft=0, hRight=0, hSize=0, hBoxMin=0, hBoxMax=0;
};

struct BenchResult {
  std::string mode;                 // "Baseline" 或 "DLB"
  KernelTimes  kt;
  std::vector<int> clusterCounts;   // 每轮簇数序列
  HashSummary  hs;
};

// 友好打印
inline void printBench(const BenchResult& r){
  printf("=== %s ===\n", r.mode.c_str());
  printf("  genNeighbours: %.3f ms\n", r.kt.genNbh);
  printf("  merge        : %.3f ms\n", r.kt.merge);
  if (r.mode=="Baseline"){
    printf("  localScan    : %.3f ms\n", r.kt.localScan);
    printf("  globalScan   : %.3f ms\n", r.kt.globalScan);
    printf("  compact(3kernels): %.3f ms\n", r.kt.compact);
  } else {
    printf("  dlbCompact(1kernel): %.3f ms\n", r.kt.compact);
  }
  printf("  TOTAL        : %.3f ms\n", r.kt.total());
  printf("  Steps        : %zu\n", r.clusterCounts.size());
  if (!r.clusterCounts.empty()){
    printf("  Start->End   : %d -> %d\n", r.clusterCounts.front(), r.clusterCounts.back());
  }
  printf("  Hash(Parents/Left/Right/Size/BoxMin/BoxMax):\n");
  printf("    %016llx %016llx %016llx %016llx %016llx %016llx\n",
    (unsigned long long)r.hs.hParents, (unsigned long long)r.hs.hLeft,
    (unsigned long long)r.hs.hRight,   (unsigned long long)r.hs.hSize,
    (unsigned long long)r.hs.hBoxMin,  (unsigned long long)r.hs.hBoxMax);
}

// 段3======= PLOC_clustering_bench.cpp (片段) =======
#include "utils_timing_hash.hpp"
#include "bench_structs.hpp"
#include "compact_dlb.cuh"

// 假定这些数组容量是 2*N - 1（BVH 全部节点）
extern int *nodeParentIndices;
extern int *nodeLeftIndices;
extern int *nodeRightIndices;
extern int *nodeSizes;
extern float3 *nodeBoxesMin;
extern float3 *nodeBoxesMax;

// 代表数组（ping-pong）
extern int *nodeIndices[2];

// 旧版三核需要：
extern int *threadOffsets;
extern int *blockOffsets;

// 生成邻居、merge 所需的参数…（你已有）

// —— 新增：DLB 资源 —— //
static int *d_blockSums = nullptr;
static int *d_blockStatus = nullptr;
static int *d_newAliveCount = nullptr;
static int   dlbCapacityBlocks = 0;

static void ensureDLBBuffers(int maxClusters){
  int needBlocks = (maxClusters + DLB_BLOCK_THREADS - 1) / DLB_BLOCK_THREADS;
  if (needBlocks <= dlbCapacityBlocks) return;
  if (d_blockSums)    cudaFree(d_blockSums);
  if (d_blockStatus)  cudaFree(d_blockStatus);
  dlbCapacityBlocks = needBlocks;
  CUDA_CHECK(cudaMalloc(&d_blockSums,   dlbCapacityBlocks * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_blockStatus, dlbCapacityBlocks * sizeof(int)));
  if (!d_newAliveCount) CUDA_CHECK(cudaMalloc(&d_newAliveCount, sizeof(int)));
}

static HashSummary makeHashSummary(int totalNodeCapacity){
  HashSummary hs{};
  // 注意：容量/有效长度按你的实现来，这里用 totalNodeCapacity 近似
  {
    auto h = d2h_copy(nodeParentIndices, totalNodeCapacity);
    hs.hParents = fnv1a64(h.data(), h.size()*sizeof(int));
  }
  {
    auto h = d2h_copy(nodeLeftIndices, totalNodeCapacity);
    hs.hLeft = fnv1a64(h.data(), h.size()*sizeof(int));
  }
  {
    auto h = d2h_copy(nodeRightIndices, totalNodeCapacity);
    hs.hRight = fnv1a64(h.data(), h.size()*sizeof(int));
  }
  {
    auto h = d2h_copy(nodeSizes, totalNodeCapacity);
    hs.hSize = fnv1a64(h.data(), h.size()*sizeof(int));
  }
  {
    auto h = d2h_copy(nodeBoxesMin, totalNodeCapacity);
    hs.hBoxMin = fnv1a64(h.data(), h.size()*sizeof(float3));
  }
  {
    auto h = d2h_copy(nodeBoxesMax, totalNodeCapacity);
    hs.hBoxMax = fnv1a64(h.data(), h.size()*sizeof(float3));
  }
  return hs;
}

// —— 包一层：单次 clustering，返回 BenchResult —— //
BenchResult runClusteringOnce(int numberOfTriangles, bool useDLB){
  BenchResult R;
  R.mode = useDLB ? "DLB" : "Baseline";
  KernelTimes &KT = R.kt;

  // 初始化
  int numberOfClusters = numberOfTriangles;
  bool swapBuffers = sortSwap; // 你的现有开关
  steps = 0;

  // DLB 资源
  if (useDLB) ensureDLBBuffers(numberOfTriangles);

  // 计时器
  CudaTimer t;

  // 主循环
  while (numberOfClusters > 1){
    ++steps;

    // === 生成邻居 ===
    t.start();
    // …… 你现有的 generateNeighboursCached / generateNeighbours 两分支 （用 launchTimed 也行）
    // （如果仍想用原作者的 launchTimed，就把返回时间加到 KT.genNbh）
    // 这里演示用事件计时包裹：
    // launchGenerateNeighbours(...);
    KT.genNbh += t.stop();

    // === 清 prefixScanOffset（可留） ===
    module->getGlobal("prefixScanOffset").clear();

    // === merge ===
    t.start();
    // launchMerge(...);   // 与原代码一致
    KT.merge += t.stop();

    // === 计算新簇数 & 压缩 ===
    if (!useDLB){
      // ---------- Baseline: 三核 ----------
      // 1) local prefix scan
      {
        int numberOfBlocks = divCeil(numberOfClusters, PLOC_SCAN_BLOCK_THREADS);
        t.start();
        // localPrefixScan<<<numberOfBlocks, PLOC_SCAN_BLOCK_THREADS>>>(...);
        KT.localScan += t.stop();
        // 2) global prefix scan
        t.start();
        // globalPrefixScan<<<1, PLOC_SCAN_BLOCK_THREADS>>>(...);
        KT.globalScan += t.stop();
        // 3) compact
        t.start();
        // compact<<<numberOfBlocks, PLOC_SCAN_BLOCK_THREADS>>>(...);
        KT.compact += t.stop();
      }
      // 新的簇数：沿用原逻辑
      int removed = 0;
      CUDA_CHECK(cudaMemcpy(&removed,
         module->getGlobal("prefixScanOffset").getPtr(), sizeof(int), cudaMemcpyDeviceToHost));
      numberOfClusters = numberOfClusters - removed;

      // 交换
      swapBuffers = !swapBuffers;
    } else {
      // ---------- DLB: 单核 ----------
      int nBlocks = divCeil(numberOfClusters, DLB_BLOCK_THREADS);
      CUDA_CHECK(cudaMemset(d_blockStatus, 0, nBlocks * sizeof(int)));
      CUDA_CHECK(cudaMemset(d_newAliveCount, 0, sizeof(int)));
      void* args[] = {
        &numberOfClusters,
        &nodeIndices[swapBuffers],
        &nodeIndices[!swapBuffers],
        &d_blockSums, &d_blockStatus, &d_newAliveCount
      };
      t.start();
      CUDA_CHECK(cudaLaunchKernel(
        (void*)compactDecoupledLookback, dim3(nBlocks), dim3(DLB_BLOCK_THREADS),
        nullptr, 0, args));
      CUDA_CHECK(cudaDeviceSynchronize());
      KT.compact += t.stop();

      int alive = 0;
      CUDA_CHECK(cudaMemcpy(&alive, d_newAliveCount, sizeof(int), cudaMemcpyDeviceToHost));
      numberOfClusters = alive;

      swapBuffers = !swapBuffers;
    }

    // 记录每一轮的簇数（用于对比）
    R.clusterCounts.push_back(numberOfClusters);
  }

  // 汇总哈希（注意容量按你的实现来；这里用 2*N-1）
  int totalNodeCapacity = 2*numberOfTriangles - 1;
  R.hs = makeHashSummary(totalNodeCapacity);
  return R;
}
// ======= main_or_driver.cpp (片段) =======
#include "bench_structs.hpp"

// 简单的一致性检查：
// 1) 步数（序列长度）相同；
// 2) 每轮的簇数一致；
// 3) 关键数组哈希一致。
static bool compareBench(const BenchResult& A, const BenchResult& B){
  bool ok = true;
  if (A.clusterCounts.size() != B.clusterCounts.size()){
    fprintf(stderr, "[DIFF] steps: %zu vs %zu\n",
      A.clusterCounts.size(), B.clusterCounts.size());
    ok = false;
  }
  size_t m = std::min(A.clusterCounts.size(), B.clusterCounts.size());
  for (size_t i=0;i<m;i++){
    if (A.clusterCounts[i] != B.clusterCounts[i]){
      fprintf(stderr, "[DIFF] round %zu clusters: %d vs %d\n",
        i, A.clusterCounts[i], B.clusterCounts[i]);
      ok = false;
      break;
    }
  }

  if (A.hs.hParents != B.hs.hParents){ fprintf(stderr,"[DIFF] parents hash\n"); ok=false; }
  if (A.hs.hLeft   != B.hs.hLeft  ){ fprintf(stderr,"[DIFF] left hash\n");    ok=false; }
  if (A.hs.hRight  != B.hs.hRight ){ fprintf(stderr,"[DIFF] right hash\n");   ok=false; }
  if (A.hs.hSize   != B.hs.hSize  ){ fprintf(stderr,"[DIFF] size hash\n");    ok=false; }
  if (A.hs.hBoxMin != B.hs.hBoxMin){ fprintf(stderr,"[DIFF] boxMin hash\n");  ok=false; }
  if (A.hs.hBoxMax != B.hs.hBoxMax){ fprintf(stderr,"[DIFF] boxMax hash\n");  ok=false; }

  return ok;
}

void runABTest(int numberOfTriangles){
  // 固定随机种子 / 输入（确保两次输入完全一致）
  // 你已有数据准备流程的话，这里只要不要在两次运行间改动输入即可

  auto baseline = runClusteringOnce(numberOfTriangles, /*useDLB=*/false);
  auto dlb      = runClusteringOnce(numberOfTriangles, /*useDLB=*/true);

  printBench(baseline);
  printBench(dlb);

  bool ok = compareBench(baseline, dlb);
  printf("=== CONSISTENCY: %s ===\n", ok ? "OK ✅" : "MISMATCH ❌");
}
