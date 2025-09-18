# main_3d.jl
# Author: Manolo Venturin
# Restart Julia to run the program

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
filename = joinpath(problem_path, "mesh2.dat")
mesh = read_mesh(filename)

filename = joinpath(problem_path, "solution2.dat")
u = read_solution(filename)

# Orient faces
# elements = orient_faces_ccw!(mesh.coordinates, mesh.elements)

# Unpack data
x = mesh.coordinates[:, 1]
y = mesh.coordinates[:, 2]
elements = mesh.elements
z = u


# Makie visualization - Mesh


# Points
points = Point3f.(x, y, 0)

# Triangles 0-based
faces = [TriangleFace(elements[i, :]) for i in 1:size(elements,1)]

# Mesh generation
mesh_data = GeometryBasics.Mesh(points, faces)

edges = Set{Tuple{Int,Int}}()
for face in faces
    i, j, k = Tuple(face)
    push!(edges, (i, j))
    push!(edges, (j, k))
    push!(edges, (k, i))
end

line_points = Point3f[]
for (i, j) in edges
    push!(line_points, points[i])
    push!(line_points, points[j])
end


# Show mesh
fig = Figure()
ax = Axis3(fig[1, 1])
mesh!(ax, mesh_data, color = :lightblue, shading = true)
linesegments!(ax, line_points, color=:black)
fig



# Makie visualization - Surf


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
