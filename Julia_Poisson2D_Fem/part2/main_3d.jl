# main_3d.jl
# Author: Manolo Venturin

using Plots
using GLMakie
using GeometryBasics

# Load mesh module
# include(joinpath(@__DIR__, "mesh.jl"))
include("mesh.jl")
using .Mesh2D

# Problem  configuration
problem_name = "Problem 4"
problem_dir = "problem4"
problem_path = joinpath(@__DIR__, "..", "data", problem_dir)

# Get data: mesh and solution
filename = joinpath(problem_path, "mesh.dat")
mesh = read_mesh(filename)

filename = joinpath(problem_path, "solution.dat")
u = read_solution(filename)

# Unpack data
x = mesh.coordinates[:, 1]
y = mesh.coordinates[:, 2]
z = u

# 3D visualization

# Points
points = Point3f.(x, y, z)

# Triangles 0-based
faces = [TriangleFace(mesh.elements[i, :]) for i in 1:size(mesh.elements,1)]

# Mesh
geom_mesh = GeometryBasics.Mesh(points, faces)

# Plot
fig = Figure()
ax = Axis3(fig[1, 1], xlabel="x", ylabel="y", zlabel="u", title="FEM solution", aspect=:data)

plt = GLMakie.mesh!(ax, geom_mesh, color = z, colormap = :viridis)

Colorbar(fig[1, 2], plt, label = "u")

fig
