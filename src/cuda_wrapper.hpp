// src/cuda_wrappers.hpp

#pragma once

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <iostream>
#include <stdexcept>

namespace cuda_wrappers {

// Error checking macro
#define CUDA_CHECK(ans) { ::cuda_wrappers::gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
        std::cerr << "GPUassert: " << cudaGetErrorString(code) << " " << file << " " << line << std::endl;
        if (abort) exit(code);
    }
}

// Memory management
inline void* malloc(size_t bytes) {
    void* ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&ptr, bytes));
    return ptr;
}
inline void free(void* ptr) {
    CUDA_CHECK(cudaFree(ptr));
}
inline void memcpy(void* dst, const void* src, size_t bytes, cudaMemcpyKind kind) {
    CUDA_CHECK(cudaMemcpy(dst, src, bytes, kind));
}

// Thrust wrappers
template<typename T>
using device_vector = thrust::device_vector<T>;
template<typename T>
using host_vector = thrust::host_vector<T>;

template<typename T>
inline void fill(device_vector<T>& vec, T value) {
    thrust::fill(vec.begin(), vec.end(), value);
}
template<typename T>
inline void copy(const device_vector<T>& src, device_vector<T>& dst) {
    thrust::copy(src.begin(), src.end(), dst.begin());
}
template<typename T>
inline void copy(const device_vector<T>& src, host_vector<T>& dst) {
    thrust::copy(src.begin(), src.end(), dst.begin());
}
template<typename T>
inline void copy(const host_vector<T>& src, device_vector<T>& dst) {
    thrust::copy(src.begin(), src.end(), dst.begin());
}

// cuSPARSE CSR matrix wrapper
struct CSRMatrix {
    int rows, cols, nnz;
    float* d_values;
    int* d_rowPtr;
    int* d_colInd;
    cusparseMatDescr_t descr;

    CSRMatrix(int r, int c, int n, float* v, int* rp, int* ci)
        : rows(r), cols(c), nnz(n), d_values(v), d_rowPtr(rp), d_colInd(ci) {
        cusparseCreateMatDescr(&descr);
        cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    }
    ~CSRMatrix() {
        cusparseDestroyMatDescr(descr);
    }
};

// Matrix-vector multiply (single precision)
inline void csr_matvec(const CSRMatrix& mat, const float* d_x, float* d_y, cusparseHandle_t handle) {
    const float alpha = 1.0f, beta = 0.0f;
    cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                   mat.rows, mat.cols, mat.nnz,
                   &alpha, mat.descr, mat.d_values, mat.d_rowPtr, mat.d_colInd,
                   d_x, &beta, d_y);
}

// cuBLAS vector operations (example: dot product)
inline float dot(int n, const float* x, int incx, const float* y, int incy, cublasHandle_t handle) {
    float result;
    cublasSdot(handle, n, x, incx, y, incy, &result);
    return result;
}

// Random number generation (stub, to be implemented as needed)
inline void init_curand_states(curandState* d_states, int n, unsigned long seed);

// CUDA streams
inline cudaStream_t create_stream() {
    cudaStream_t stream;
    CUDA_CHECK(cudaStreamCreate(&stream));
    return stream;
}
inline void destroy_stream(cudaStream_t stream) {
    CUDA_CHECK(cudaStreamDestroy(stream));
}
inline void stream_synchronize(cudaStream_t stream) {
    CUDA_CHECK(cudaStreamSynchronize(stream));
}

// CUDA events
inline cudaEvent_t create_event() {
    cudaEvent_t event;
    CUDA_CHECK(cudaEventCreate(&event));
    return event;
}
inline void destroy_event(cudaEvent_t event) {
    CUDA_CHECK(cudaEventDestroy(event));
}
inline float elapsed_time(cudaEvent_t start, cudaEvent_t end) {
    float ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&ms, start, end));
    return ms;
}

} // namespace cuda_wrappers