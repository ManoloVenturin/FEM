# main_mtx_operations.jl

using SparseArrays
using MatrixDepot
using Plots

# Create a Poisson matrix of size (n^2 x n^2)
n = 10
A = matrixdepot("poisson", n)

# Print matrix info
println("Matrix size: ", size(A))
println("NNZ: ", nnz(A))

# Show matrix pattern
spy(A, title="Sparsity Pattern - Poisson Matrix (n = $n)", markersize=1)
