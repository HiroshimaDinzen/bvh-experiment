// morton_sort_benchmark.cu
// nvcc -O3 -std=c++17 morton_sort_benchmark.cu -o morton_bench
// 运行示例：
//   ./morton_bench               // 默认 N=4M, BINS=1024, streams=8
//   ./morton_bench 8000000 2048  // N=8M, BINS=2048, streams=8
//   ./morton_bench 2000000 512 4 // N=2M, BINS=512, streams=4
//
// 说明：
//   方法A：全局一次 radix sort (Thrust/CUB)
//   方法B：粗分桶(按高位) + 桶内 radix sort（分桶稳定线性，桶内排序可并行）
//   校验：方法B的结果与方法A逐元素相同视为通过
//   输出：两种方法的耗时(ms)与加速比
//
// 注意：请确保已安装 CUDA 11+（建议 12/13），并可用 Thrust/CUB。

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <iostream>

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

// ------------------------ CUDA error check ------------------------
#define CUDA_CHECK(expr) \
do { \
    cudaError_t __err = (expr); \
    if (__err != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d : %s\n", #expr, __FILE__, __LINE__, cudaGetErrorString(__err)); \
        std::exit(1); \
    } \
} while(0)

// ------------------------ Morton encode (30-bit) ------------------
// 扩展10位到30位（插零）参考经典 bit interleave 实现
__host__ __device__ inline uint32_t part1by2(uint32_t x) {
    x &= 0x000003ffu;                 // 10 bits
    x = (x ^ (x << 16)) & 0x030000FFu;
    x = (x ^ (x << 8 )) & 0x0300F00Fu;
    x = (x ^ (x << 4 )) & 0x030C30C3u;
    x = (x ^ (x << 2 )) & 0x09249249u;
    return x;
}

// 把每轴 10bit 整数编码为 30-bit Morton（XYZ 交织）
__host__ __device__ inline uint32_t morton3D_10bits(uint32_t xi, uint32_t yi, uint32_t zi) {
    return (part1by2(zi) << 2) | (part1by2(yi) << 1) | (part1by2(xi));
}

// 将 [min,max] 归一化到 [0,1] 并量化为 10-bit 整数（0..1023）
__host__ __device__ inline uint32_t quant10(float v, float vmin, float vmax) {
    float t = (v - vmin) / (vmax - vmin);
    t = fminf(fmaxf(t, 0.0f), 1.0f);
    const float scale = 1023.0f; // 2^10 - 1
    return static_cast<uint32_t>(t * scale + 0.5f);
}

// ------------------------ Kernels ------------------------

// 由 Morton key 计算桶号：取高位 topBits（BINS=2^topBits）
__global__ void compute_bin_ids(const uint32_t* __restrict__ keys,
                                uint32_t* __restrict__ bin_ids,
                                uint32_t N, int topBits)
{
    uint32_t i = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (i >= N) return;
    uint32_t k = keys[i];
    bin_ids[i] = (topBits == 0) ? 0u : (k >> (32 - topBits));
}

// 原地直方图（原子加）
__global__ void histogram_bins(const uint32_t* __restrict__ bin_ids,
                               uint32_t* __restrict__ hist, uint32_t N)
{
    uint32_t i = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (i >= N) return;
    atomicAdd(hist + bin_ids[i], 1u);
}

// 按桶 Scatter：根据偏移表把 (key,val) 写到分桶后的有序区间
__global__ void scatter_by_bins(const uint32_t* __restrict__ bin_ids,
                                const uint32_t* __restrict__ in_keys,
                                const uint32_t* __restrict__ in_vals,
                                uint32_t* __restrict__ out_keys,
                                uint32_t* __restrict__ out_vals,
                                uint32_t* __restrict__ bin_write_ptrs,
                                uint32_t N)
{
    uint32_t i = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (i >= N) return;
    uint32_t b = bin_ids[i];
    // 原子增加，得到当前元素在该桶的全局写入位置
    uint32_t pos = atomicAdd(bin_write_ptrs + b, 1u);
    out_keys[pos] = in_keys[i];
    out_vals[pos] = in_vals[i];
}

// ------------------------ Utils ------------------------
static inline bool is_power_of_two(uint32_t x) {
    return x && ((x & (x - 1)) == 0);
}

static int ilog2_uint(uint32_t x) {
    int r = -1;
    while (x) { x >>= 1; ++r; }
    return r;
}

// 生成 N 个随机点并计算 Morton key（host 端生成）
void generate_morton_keys(thrust::host_vector<uint32_t>& hkeys,
                          thrust::host_vector<uint32_t>& hvals,
                          uint32_t N,
                          uint64_t seed = 20250917ULL)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<float> uni(0.0f, 1.0f);

    // 随机生成三维点并映射到 [xmin,xmax] 区间（可调）
    const float xmin = -100.0f, xmax = 300.0f;
    const float ymin = -50.0f,  ymax = 200.0f;
    const float zmin = -10.0f,  zmax = 510.0f;

    hkeys.resize(N);
    hvals.resize(N);

    for (uint32_t i = 0; i < N; ++i) {
        float x = xmin + (xmax - xmin) * uni(rng);
        float y = ymin + (ymax - ymin) * uni(rng);
        float z = zmin + (zmax - zmin) * uni(rng);
        uint32_t xi = quant10(x, xmin, xmax);
        uint32_t yi = quant10(y, ymin, ymax);
        uint32_t zi = quant10(z, zmin, zmax);
        uint32_t key = morton3D_10bits(xi, yi, zi);
        hkeys[i] = key;
        hvals[i] = i; // 记录原索引
    }
}

// 校验两份(key,val)是否完全一致
bool check_equal(const thrust::host_vector<uint32_t>& a_keys,
                 const thrust::host_vector<uint32_t>& a_vals,
                 const thrust::host_vector<uint32_t>& b_keys,
                 const thrust::host_vector<uint32_t>& b_vals)
{
    if (a_keys.size() != b_keys.size()) return false;
    if (a_vals.size() != b_vals.size()) return false;
    size_t N = a_keys.size();
    for (size_t i = 0; i < N; ++i) {
        if (a_keys[i] != b_keys[i] || a_vals[i] != b_vals[i]) {
            fprintf(stderr, "Mismatch at %zu: (%u,%u) vs (%u,%u)\n",
                    i, a_keys[i], a_vals[i], b_keys[i], b_vals[i]);
            return false;
        }
    }
    return true;
}

// ------------------------ Main Benchmark ------------------------
int main(int argc, char** argv)
{
    // 参数
    uint32_t N = 4u * 1024u * 1024u; // 默认 4M
    uint32_t BINS = 1024;            // 默认 1024 桶（取 key 的高 10 bit）
    int num_streams = 8;             // 默认 8 个 stream 并行桶内排序

    if (argc >= 2) N = static_cast<uint32_t>(std::stoul(argv[1]));
    if (argc >= 3) BINS = static_cast<uint32_t>(std::stoul(argv[2]));
    if (argc >= 4) num_streams = std::max(1, std::stoi(argv[3]));

    if (!is_power_of_two(BINS) || BINS > (1u << 16)) {
        fprintf(stderr, "BINS 必须为 2 的幂且不超过 65536，当前 BINS=%u\n", BINS);
        return 1;
    }
    int topBits = ilog2_uint(BINS);
    printf("Config: N=%u, BINS=%u (topBits=%d), streams=%d\n", N, BINS, topBits, num_streams);

    // 生成随机 Morton key
    thrust::host_vector<uint32_t> hkeys, hvals;
    generate_morton_keys(hkeys, hvals, N);

    // 拷到 GPU
    thrust::device_vector<uint32_t> d_keys = hkeys;
    thrust::device_vector<uint32_t> d_vals = hvals;

    // -------- 方法A：全局一次 radix sort --------
    thrust::device_vector<uint32_t> A_keys = d_keys;
    thrust::device_vector<uint32_t> A_vals = d_vals;

    cudaEvent_t eA_start, eA_stop;
    CUDA_CHECK(cudaEventCreate(&eA_start));
    CUDA_CHECK(cudaEventCreate(&eA_stop));
    CUDA_CHECK(cudaEventRecord(eA_start));

    thrust::sort_by_key(A_keys.begin(), A_keys.end(), A_vals.begin()); // 全局一次

    CUDA_CHECK(cudaEventRecord(eA_stop));
    CUDA_CHECK(cudaEventSynchronize(eA_stop));
    float tA_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&tA_ms, eA_start, eA_stop));

    // -------- 方法B：粗分桶 + 桶内排序 --------
    thrust::device_vector<uint32_t> B_in_keys = d_keys;
    thrust::device_vector<uint32_t> B_in_vals = d_vals;
    thrust::device_vector<uint32_t> B_bin_ids(N);
    thrust::device_vector<uint32_t> B_hist(BINS, 0u);
    thrust::device_vector<uint32_t> B_offsets(BINS, 0u);
    thrust::device_vector<uint32_t> B_write_ptrs(BINS, 0u);
    thrust::device_vector<uint32_t> B_grouped_keys(N);
    thrust::device_vector<uint32_t> B_grouped_vals(N);

    const int block = 256;
    const int gridN = (N + block - 1) / block;

    cudaEvent_t eB_start, eB_mid, eB_stop;
    CUDA_CHECK(cudaEventCreate(&eB_start));
    CUDA_CHECK(cudaEventCreate(&eB_mid));
    CUDA_CHECK(cudaEventCreate(&eB_stop));

    CUDA_CHECK(cudaEventRecord(eB_start));

    // 1) 计算桶号
    compute_bin_ids<<<gridN, block>>>(thrust::raw_pointer_cast(B_in_keys.data()),
                                      thrust::raw_pointer_cast(B_bin_ids.data()),
                                      N, topBits);
    CUDA_CHECK(cudaGetLastError());

    // 2) 直方图
    histogram_bins<<<gridN, block>>>(thrust::raw_pointer_cast(B_bin_ids.data()),
                                     thrust::raw_pointer_cast(B_hist.data()), N);
    CUDA_CHECK(cudaGetLastError());

    // 3) 前缀和 => 每个桶的起始偏移
    thrust::exclusive_scan(B_hist.begin(), B_hist.end(), B_offsets.begin());

    // 4) 初始化写指针 = offsets
    thrust::copy(B_offsets.begin(), B_offsets.end(), B_write_ptrs.begin());

    // 5) Scatter 到分桶后的连续区间
    scatter_by_bins<<<gridN, block>>>(
        thrust::raw_pointer_cast(B_bin_ids.data()),
        thrust::raw_pointer_cast(B_in_keys.data()),
        thrust::raw_pointer_cast(B_in_vals.data()),
        thrust::raw_pointer_cast(B_grouped_keys.data()),
        thrust::raw_pointer_cast(B_grouped_vals.data()),
        thrust::raw_pointer_cast(B_write_ptrs.data()),
        N
    );
    CUDA_CHECK(cudaGetLastError());

    // 记录到“分桶完毕”时刻
    CUDA_CHECK(cudaEventRecord(eB_mid));
    CUDA_CHECK(cudaEventSynchronize(eB_mid));

    // 6) 桶内排序（每个 [offset[i], offset[i+1]) 区间 sort_by_key）
    //    可使用多个 stream 提高并行度
    std::vector<cudaStream_t> streams(num_streams);
    for (int s = 0; s < num_streams; ++s) CUDA_CHECK(cudaStreamCreate(&streams[s]));

    // 拉下 offsets 到 host 以便划分区间
    thrust::host_vector<uint32_t> h_offsets = B_offsets;
    thrust::host_vector<uint32_t> h_counts(BINS);
    for (uint32_t i = 0; i < BINS; ++i) {
        uint32_t start = h_offsets[i];
        uint32_t end   = (i + 1 < BINS) ? h_offsets[i + 1] : N;
        h_counts[i] = end - start;
    }

    // 调度：按桶大小降序可更好利用并行，这里简单按序分配
    for (uint32_t b = 0; b < BINS; ++b) {
        uint32_t start = h_offsets[b];
        uint32_t cnt   = h_counts[b];
        if (cnt <= 1) continue;
        int sid = (num_streams > 0) ? (b % num_streams) : 0;
        auto k_begin = B_grouped_keys.begin() + start;
        auto k_end   = k_begin + cnt;
        auto v_begin = B_grouped_vals.begin() + start;
        // 使用带流的 thrust execution policy
        thrust::sort_by_key(thrust::cuda::par.on(streams[sid]), k_begin, k_end, v_begin);
    }
    // 等待所有流完成
    for (int s = 0; s < num_streams; ++s) CUDA_CHECK(cudaStreamSynchronize(streams[s]));

    CUDA_CHECK(cudaEventRecord(eB_stop));
    CUDA_CHECK(cudaEventSynchronize(eB_stop));
    float tB_total_ms = 0.0f, tB_bucket_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&tB_total_ms, eB_start, eB_stop));
    CUDA_CHECK(cudaEventElapsedTime(&tB_bucket_ms, eB_start, eB_mid));

    // 清理流
    for (int s = 0; s < num_streams; ++s) CUDA_CHECK(cudaStreamDestroy(streams[s]));

    // 方法B的分桶+桶内排序结果即为最终全局有序
    thrust::host_vector<uint32_t> hA_keys = A_keys;
    thrust::host_vector<uint32_t> hA_vals = A_vals;
    thrust::host_vector<uint32_t> hB_keys = B_grouped_keys;
    thrust::host_vector<uint32_t> hB_vals = B_grouped_vals;

    bool ok = check_equal(hA_keys, hA_vals, hB_keys, hB_vals);

    // 输出结果
    printf("\n===== Benchmark Result =====\n");
    printf("Method A (Global radix sort): %.3f ms\n", tA_ms);
    printf("Method B (Coarse binning + per-bin sort): total=%.3f ms (binning=%.3f ms, per-bin-sort=%.3f ms)\n",
           tB_total_ms, tB_bucket_ms, tB_total_ms - tB_bucket_ms);
    if (tB_total_ms > 0.0f)
        printf("Speedup (A/B): %.3f x\n", tA_ms / tB_total_ms);
    printf("Correctness: %s\n", ok ? "PASS" : "FAIL");

    // 资源回收
    CUDA_CHECK(cudaEventDestroy(eA_start));
    CUDA_CHECK(cudaEventDestroy(eA_stop));
    CUDA_CHECK(cudaEventDestroy(eB_start));
    CUDA_CHECK(cudaEventDestroy(eB_mid));
    CUDA_CHECK(cudaEventDestroy(eB_stop));

    return ok ? 0 : 2;
}
