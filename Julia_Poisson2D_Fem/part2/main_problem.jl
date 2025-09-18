# main_problem.jl
# Author: Manolo Venturin
# Restart Julia to run the program

using Plots

# Load mesh module
# include(joinpath(@__DIR__, "mesh.jl"))
include("mesh.jl")
using .Mesh2D

# Problem configuration
# problem_name = "Problem 1"
# problem_name = "Problem 2"
# problem_name = "Problem 3"
problem_name = "Problem 4"

if problem_name == "Problem 1"
    problem_dir = "problem1"
    bcstyles = Dict(
        1 => (:red, :solid),
        2 => (:red, :solid),
        3 => (:red, :solid),
        4 => (:red, :solid)
    )
end

if problem_name == "Problem 2"
    problem_dir = "problem2"
    bcstyles = Dict(
        1 => (:red, :solid),
        2 => (:red, :solid),
        3 => (:red, :solid),
        4 => (:red, :solid)
    )
end

if problem_name == "Problem 3"
    problem_dir = "problem3"
    bcstyles = Dict(
        1 => (:red, :solid),
        2 => (:blue, :solid)
    )
end

if problem_name == "Problem 4"
    problem_dir = "problem4"
    bcstyles = Dict(
        1 => (:red, :solid),
        2 => (:blue, :solid),
        3 => (:green, :solid)
    )
end

problem_path = joinpath(@__DIR__, "..", "data", problem_dir)
output_dir = joinpath(@__DIR__, "img")

# Problem
println("=== $(problem_name) ===\n")

# Get data: mesh and solution
filename = joinpath(problem_path, "mesh.dat")
mesh = read_mesh(filename)

filename = joinpath(problem_path, "solution.dat")
u = read_solution(filename)

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
plot_solution(mesh, u)

filename = joinpath(@__DIR__, "img", problem_dir * "_solution1.png")
savefig(filename)

# This code has a bug in tricontour!
# # Plot (contour)
# plot_solution(mesh, u; style="contour")

# filename = joinpath(@__DIR__, "img", problem_dir * "_solution2.png")
# savefig(filename)
