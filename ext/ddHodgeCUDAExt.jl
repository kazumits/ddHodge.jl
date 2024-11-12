module ddHodgeCUDAExt

using ddHodge.GPU
using CUDA: CuVector, CUSPARSE.CuSparseMatrixCSR

GPU.cuspm(X) = CuSparseMatrixCSR(X)
GPU.cuvec(x) = CuVector(x)

end # module ddHodgeCUDAExt

