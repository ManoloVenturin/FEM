using LinearAlgebra
using MatrixDepot
using SparseArrays

function steepest_descent(A, b; x0=zeros(length(b)), tol=1e-8, maxiter=1000)
    x = copy(x0)
    r = b - A * x
    b_norm = norm(b)
    rel_residuals = [norm(r) / b_norm]

    for k in 1:maxiter
        if norm(r) < tol
            break
        end
        alpha = dot(r, r) / dot(r, A * r)
        x += alpha * r
        r = b - A * x
        push!(rel_residuals, norm(r) / b_norm)
    end
    return x, rel_residuals
end

# Create Poisson matrix of size (n^2 × n^2)
n = 10
A = matrixdepot("poisson", n)

# Create a reference solution and vector b
x_true = ones(size(A, 2))
b = A * x_true

# Print matrix info
println("Matrix size: ", size(A))
println("NNZ: ", nnz(A))

# Main
x, residuals = steepest_descent(A, b)
plot(residuals, yscale = :log10, xlabel = "Iteration", ylabel = "Relative residual", title = "SD", legend = false)
