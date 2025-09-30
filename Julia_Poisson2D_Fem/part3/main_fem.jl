# main_fem.jl
# Author: Manolo Venturin
# Restart Julia to run the program

using DelimitedFiles
using Plots

# Load mesh module
# include(joinpath(@__DIR__, "mesh.jl"))
include("mesh.jl")
using .Mesh2D

# Problem configuration
problem_name = "Problem 4"
problem_dir = "problem4"
bcstyles = Dict(
    1 => (:red, :solid),
    2 => (:blue, :solid),
    3 => (:green, :solid)
)

problem_path = joinpath(@__DIR__, "..", "data", problem_dir)
output_dir = joinpath(@__DIR__, "img")

# Problem
println("=== $(problem_name) ===\n")

# Get data: mesh and solution
filename = joinpath(problem_path, "mesh.dat")
mesh = read_mesh(filename)

filename = joinpath(problem_path, "solution.dat")
u_ref = read_solution(filename)

# Matrix data

function read_freefem_vector(filename::String)
    # Apri il file
    lines = readlines(filename)

    # Leggi il numero di valori dalla prima riga
    n_values = parse(Int, strip(lines[1]))

    # Leggi tutte le altre righe e dividile in parole (numeri)
    values = Float64[]

    for line in lines[2:end]
        # Rimuove spazi multipli e dividi in parole
        tokens = split(strip(line))

        # Converti ogni parola in Float64 e aggiungila alla lista
        append!(values, parse.(Float64, tokens))
    end

    # Verifica che abbiamo letto il numero giusto di valori
    if length(values) != n_values
        error("Expected $n_values values, but read $(length(values))")
    end

    return values
end

filename = joinpath(problem_path, "F_data.dat")
F_data = read_freefem_vector(filename)

# Basic info
print_mesh_info(mesh)

# Show mesh data for debug
if false
    @show mesh.coordinates
    @show mesh.elements
    @show mesh.boundary_edges
    @show mesh.boundary_tags
end

# Domain
plot_domain(mesh; tag_styles=bcstyles, showtag=true)

filename = joinpath(output_dir, problem_dir * "_domain.png")
savefig(filename)

# Mesh
plot_mesh(mesh)

filename = joinpath(output_dir, problem_dir * "_mesh.png")
savefig(filename)

# Plot (fill)
plot_solution(mesh, u_ref)

filename = joinpath(@__DIR__, "img", problem_dir * "_solution1.png")
savefig(filename)

# FEM
module MyPoisson

function get_a(Omega_e, x, y)
    i, j, k = Omega_e
    a_i = x[j] * y[k] - x[k] * y[j]
    a_j = x[k] * y[i] - x[i] * y[k]
    a_k = x[i] * y[j] - x[j] * y[i]
    return (a_i, a_j, a_k)
end

function get_b(Omega_e, x, y)
    i, j, k = Omega_e
    b_i = y[j] - y[k]
    b_j = y[k] - y[i]
    b_k = y[i] - y[j]
    return (b_i, b_j, b_k)
end

function get_c(Omega_e, x, y)
    i, j, k = Omega_e
    c_i = x[k] - x[j]
    c_j = x[i] - x[k]
    c_k = x[j] - x[i]
    return (c_i, c_j, c_k)
end

function get_A(Omega_e, x, y)
    i, j, k = Omega_e
    A2 = (x[i] * y[j] - x[j] * y[i]) + (x[k] * y[i] - x[i] * y[k]) + (x[j] * y[k] - x[k] * y[j])
    return abs(A2) / 2.0
end

function get_total_area(mesh)
    coordinates = mesh.coordinates
    x = coordinates[:,1]
    y = coordinates[:,2]

    TotalArea = 0.
    for ind_e=1:size(mesh.elements,1)
        Omega_e = mesh.elements[ind_e,:] 
        A = get_A(Omega_e, x, y)
        TotalArea += A
    end
    return TotalArea
end

function get_Ke(Omega_e, x, y)
    b = get_b(Omega_e, x, y)
    c = get_c(Omega_e, x, y)
    A = get_A(Omega_e, x, y)

    Ke = zeros(3,3)
    # Ke = (1/(4*A)) * ([b_i, b_j, b_k]' * [b_i, b_j, b_k] + [c_i, c_j, c_k]' * [c_i, c_j, c_k])
    for i=1:3
        for j=1:3
            Ke[i,j] = (b[i] * b[j] + c[i] * c[j]) / (4.0*A)
        end
    end
    return Ke
end

function assemble_K(mesh)
    coordinates = mesh.coordinates
    x = coordinates[:,1]
    y = coordinates[:,2]

    n = size(mesh.coordinates,1)
    K = zeros(n, n)
    for ind_e=1:size(mesh.elements,1)
        Omega_e = mesh.elements[ind_e,:] 
        Ke_loc = get_Ke(Omega_e, x, y)
        for i=1:3
            for j=1:3
                K[Omega_e[i],Omega_e[j]] += Ke_loc[i,j]
            end
        end        
    end
    return K
end

function get_Fe(Omega_e, x, y)
    A = get_A(Omega_e, x, y)

    fe = zeros(3)
    # fe = (A/3) * ones(3)
    for i=1:3
        fe[i] = A/3
    end
    return fe
end

function get_neumann(edge, x, y; h=1.0)
    i, j = edge
    # L = hypot(x[j]-x[i], y[j]-y[i])
    L = sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2)
    fe = (h * L / 2) * ones(2)
    return fe
end

function assemble_F(mesh)
    coordinates = mesh.coordinates
    x = coordinates[:,1]
    y = coordinates[:,2]

    n = size(mesh.coordinates,1)

    F = zeros(n)

    # Internal elements
    for ind_e=1:size(mesh.elements,1)
        Omega_e = mesh.elements[ind_e,:] 
        Fe_loc = get_Fe(Omega_e, x, y)
        for i=1:3
            F[Omega_e[i]] += Fe_loc[i]
        end        
    end

    # Boundary edges - Newumann for bc tag = 2
    for (index, tag) in enumerate(mesh.boundary_tags)
        if tag == 2
            edge = mesh.boundary_edges[index,:]
            fe_loc = get_neumann(edge, x, y; h=1.0)
            F[edge[1]] += fe_loc[1]
            F[edge[2]] += fe_loc[2]
        end
    end

    return F
end

function apply_dirichlet!(K, F, dirichlet_nodes, uD)
    for (i, value) in zip(dirichlet_nodes, uD)
        K[i, :] .= 0.0          # azzera riga
        K[:, i] .= 0.0          # azzera colonna
        K[i,i] = 1.0            # diagonale = 1
        F[i] = value            # termine noto
    end
end

function solve(mesh)
    # Dirichlet nodes with their values
    dirichlet_nodes = Int[]
    for (index, tag) in enumerate(mesh.boundary_tags)
        if tag == 1
            edge = mesh.boundary_edges[index,:]
            append!(dirichlet_nodes, edge...) 
        end
    end
    dirichlet_nodes = unique(dirichlet_nodes)
    uD_values = zeros(length(dirichlet_nodes))  # or function

    # Compute RHS and LHS
    K = assemble_K(mesh)
    F = assemble_F(mesh)

    # Apply Dirichlet boundary conditions
    apply_dirichlet!(K, F, dirichlet_nodes, uD_values)

    # Solve the linear system
    u = K \ F

    return K, F, u
end

end


TotalArea = MyPoisson.get_total_area(mesh)
println("Total area: {$TotalArea} (compare to 1*1 - π *0.3^2 ≈ 0.7173)")

# Solve the problem
K, F, u = MyPoisson.solve(mesh)


# plot_solution(mesh, u_ref)
plot_solution(mesh, u)
# plot_solution(mesh, u - u_ref)

# Compare the solution with a reference solution
u_diff = abs.(u - u_ref)
println(maximum(u_diff))

# plot_solution(mesh, u - u_ref)
