# main_mtx_linearsystem.jl

using LinearAlgebra
using MatrixDepot
using SparseArrays
using Plots
using IterativeSolvers
using IncompleteLU
using AlgebraicMultigrid

# Create Poisson matrix of size (n^2 × n^2)
n = 10
A = matrixdepot("poisson", n)

# Create a reference solution and vector b
x_true = ones(size(A, 2))
b = A * x_true

# Print matrix info
println("Matrix size: ", size(A))
println("NNZ: ", nnz(A))

# Show matrix pattern
spy(A, title="Sparsity Pattern - Poisson Matrix (n = $n)", markersize=1)

# Solving using conjugate gradient (CG)
x, history1 = cg(A, b; reltol = 1e-8, maxiter = 1000, log = true)

println("Convergence: ", history1.isconverged)
println("Iterations: ", history1.iters)
println("Final residual: ", history1.data[:resnorm][end])
residuals = history1.data[:resnorm]
plot(residuals, yscale = :log10, xlabel = "Iteration", ylabel = "Relative residual", title = "CG", legend = false)

# Jacobi precoditioner (diagonal)
M = Diagonal(1.0 ./ Vector(diag(A))) 
x, history2 = cg(A, b; Pl=M, reltol = 1e-8, maxiter = 1000, log = true)

println("Convergence: ", history2.isconverged)
println("Iterations: ", history2.iters)
println("Final residual: ", history2.data[:resnorm][end])
residuals = history2.data[:resnorm]
plot(residuals, yscale = :log10, xlabel = "Iteration", ylabel = "Relative residual", title = "CG", legend = false)

# Incomplete Cholesky
ic = ilu(A; τ=1e-1)
x, history3 = cg(A, b; Pl=ic, reltol = 1e-8, maxiter = 1000, log = true)

println("Convergence: ", history3.isconverged)
println("Iterations: ", history3.iters)
println("Final residual: ", history3.data[:resnorm][end])
residuals = history3.data[:resnorm]
plot(residuals, yscale = :log10, xlabel = "Iteration", ylabel = "Relative residual", title = "CG-IC", legend = false)

# AMG
ml = ruge_stuben(A)
p = aspreconditioner(ml)
x, history4 = cg(A, b; Pl=p, reltol = 1e-8, maxiter = 1000, log = true)

println("Convergence: ", history4.isconverged)
println("Iterations: ", history4.iters)
println("Final residual: ", history4.data[:resnorm][end])
residuals = history4.data[:resnorm]
plot(residuals, yscale = :log10, xlabel = "Iteration", ylabel = "Relative residual", title = "CG-AMG", legend = false)
