module GPUExt

using ddHodge.GPU, CUDA
GPU.cuspm(X) = CUSPARSE.CuSparseMatrixCSR(X)
GPU.cuvec(x) = CuVector(x)

end # module GPUExt
