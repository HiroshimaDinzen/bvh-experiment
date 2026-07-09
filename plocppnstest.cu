#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cuda_runtime.h>

#define SEARCH_RADIUS 4
#define ITERS 1  // 大数据可以减少循环次数避免过久

#define CUDA_CHECK(call) do { \
    cudaError_t err__ = (call); \
    if (err__ != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(err__), __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

// 欧氏距离平方（只用x有效，y,z=0）
__device__ float distance_sq(const float3& a, const float3& b) {
    float dx = a.x - b.x;
    return dx * dx; // y,z为0，不参与
}

// 单向：只向右
__global__ void nn_right_kernel(const float3* __restrict__ clusters,
                                int* __restrict__ nearestIdx, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    float minDist = 1e30f;
    int minIdx = -1;

    #pragma unroll
    for (int r = 1; r <= SEARCH_RADIUS; ++r) {
        int j = i + r;
        if (j >= n) break;
        float d = distance_sq(clusters[i], clusters[j]);
        if (d < minDist) { minDist = d; minIdx = j; }
    }
    nearestIdx[i] = minIdx;
}

// 双向：左右各 SEARCH_RADIUS
__global__ void nn_bidir_kernel(const float3* __restrict__ clusters,
                                int* __restrict__ nearestIdx, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    float minDist = 1e30f;
    int minIdx = -1;

    #pragma unroll
    for (int r = 1; r <= SEARCH_RADIUS; ++r) {
        int j1 = i + r;
        if (j1 < n) {
            float d = distance_sq(clusters[i], clusters[j1]);
            if (d < minDist) { minDist = d; minIdx = j1; }
        }
        int j2 = i - r;
        if (j2 >= 0) {
            float d = distance_sq(clusters[i], clusters[j2]);
            if (d < minDist) { minDist = d; minIdx = j2; }
        }
    }
    nearestIdx[i] = minIdx;
}

int main() {
    const int n = 1 << 20; // 2^20 = 1048576
    const float RANGE = (float)(1UL << 30); // [0, 2^30)
    float3* h_clusters = (float3*)malloc(n * sizeof(float3));
    int* h_idx_right = (int*)malloc(n * sizeof(int));
    int* h_idx_bidir = (int*)malloc(n * sizeof(int));

    // 随机生成 1D 点
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; ++i) {
        h_clusters[i].x = (float)rand() / RAND_MAX * RANGE;
        h_clusters[i].y = 0.0f;
        h_clusters[i].z = 0.0f;
    }

    // 设备内存
    float3 *d_clusters = nullptr;
    int *d_idx_right = nullptr, *d_idx_bidir = nullptr;
    CUDA_CHECK(cudaMalloc(&d_clusters, n * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&d_idx_right, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_idx_bidir, n * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_clusters, h_clusters, n * sizeof(float3), cudaMemcpyHostToDevice));

    int block = 256, grid = (n + block - 1) / block;

    // CUDA events
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    // Warm-up
    nn_right_kernel<<<grid, block>>>(d_clusters, d_idx_right, n);
    nn_bidir_kernel<<<grid, block>>>(d_clusters, d_idx_bidir, n);
    CUDA_CHECK(cudaDeviceSynchronize());

    // 计时：单向
    CUDA_CHECK(cudaEventRecord(start));
    for (int it = 0; it < ITERS; ++it) {
        nn_right_kernel<<<grid, block>>>(d_clusters, d_idx_right, n);
    }
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float ms_right_total = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&ms_right_total, start, stop));
    float ms_right_avg = ms_right_total / ITERS;

    // 计时：双向
    CUDA_CHECK(cudaEventRecord(start));
    for (int it = 0; it < ITERS; ++it) {
        nn_bidir_kernel<<<grid, block>>>(d_clusters, d_idx_bidir, n);
    }
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float ms_bidir_total = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&ms_bidir_total, start, stop));
    float ms_bidir_avg = ms_bidir_total / ITERS;

    // 拷回部分结果检查
    CUDA_CHECK(cudaMemcpy(h_idx_right, d_idx_right, n * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_idx_bidir, d_idx_bidir, n * sizeof(int), cudaMemcpyDeviceToHost));

    printf("== Indices preview (first 10) ==\n");
    printf("i\tright\tbidir\n");
    for (int i = 0; i < 10; ++i) {
        printf("%d\t%d\t%d\n", i, h_idx_right[i], h_idx_bidir[i]);
    }

    printf("\n== Timing (avg over %d iters) ==\n", ITERS);
    printf("Right-only   : %.6f ms\n", ms_right_avg);
    printf("Bidirectional: %.6f ms\n", ms_bidir_avg);
    if (ms_bidir_avg > 0.0f) {
        printf("Speedup (right / bidir): %.3fx\n", ms_right_avg / ms_bidir_avg);
    }

    // 清理
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    CUDA_CHECK(cudaFree(d_clusters));
    CUDA_CHECK(cudaFree(d_idx_right));
    CUDA_CHECK(cudaFree(d_idx_bidir));
    free(h_clusters);
    free(h_idx_right);
    free(h_idx_bidir);
    return 0;
}
