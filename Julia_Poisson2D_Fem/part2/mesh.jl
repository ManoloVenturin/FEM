# mesh.jl
# Author: Manolo Venturin

module Mesh2D

export Mesh, read_mesh, read_solution, print_mesh_info
export plot_domain, plot_mesh, plot_solution

using DelimitedFiles
using LinearAlgebra
using Plots, TriplotRecipes

"""
Data structure for a 2D triangular mesh.
"""
struct Mesh
    coordinates::Matrix{Float64}    # Node coordinates: N x 2
    elements::Matrix{Int}           # Triangle connectivity: T x 3
    boundary_edges::Matrix{Int}     # Boundary edges: B x 2
    boundary_tags::Vector{Int}      # Boundary edge labels: B
end

"""
Reads a 2D triangular mesh from a plain-text file.

# Expected file format:

The mesh file must follow this structure:

<dimension> <num_nodes> <num_elements> <num_boundary_edges>
Node coordinates
x₁ y₁
x₂ y₂
...
Triangle elements (0-based indices)
i₁₁ i₁₂ i₁₃
...
Boundary edges and tags (0-based indices)
j₁ j₂ tag
...

Note: All indices in the file are **0-based** and will be converted to **1-based** internally.
"""
function read_mesh(filename::String)::Mesh
    open(filename, "r") do io
        # --- Read header line: dimension, counts
        header_fields = split(readline(io))
        dimension             = parse(Int, header_fields[1])
        num_nodes             = parse(Int, header_fields[2])
        num_elements          = parse(Int, header_fields[3])
        num_boundary_edges    = parse(Int, header_fields[4])

        # Only 2D meshes are supported
        if dimension != 2
            error("Only 2D meshes are supported. Found dimension = $dimension")
        end

        # --- Read node coordinates
        readline(io)  # skip comment or blank line
        coordinates = zeros(Float64, num_nodes, 2)
        for node_index in 1:num_nodes
            node_fields = split(strip(readline(io)))
            coordinates[node_index, 1] = parse(Float64, node_fields[1])
            coordinates[node_index, 2] = parse(Float64, node_fields[2])
        end

        # --- Read triangle elements
        readline(io)  # skip comment or blank line
        elements = zeros(Int, num_elements, 3)
        for elem_index in 1:num_elements
            vertex_indices = parse.(Int, split(strip(readline(io))))
            elements[elem_index, :] = vertex_indices .+ 1  # 0-based -> 1-based
        end

        # --- Read boundary edges
        readline(io)  # skip comment or blank line
        boundary_edges = zeros(Int, num_boundary_edges, 2)
        boundary_tags  = zeros(Int, num_boundary_edges)
        for edge_index in 1:num_boundary_edges
            edge_fields = split(strip(readline(io)))
            boundary_edges[edge_index, :] = parse.(Int, edge_fields[1:2]) .+ 1
            boundary_tags[edge_index]     = parse(Int, edge_fields[3])
        end

        return Mesh(coordinates, elements, boundary_edges, boundary_tags)
    end
end

"""
Read nodal scalar solution u from a text file.
"""
function read_solution(filename::String)
    data = readdlm(filename)
    return data[:, 3]
end

"""
Prints a summary of the mesh.
"""
function print_mesh_info(mesh::Mesh)
    num_nodes = size(mesh.coordinates, 1)
    num_elements = size(mesh.elements, 1)
    num_boundary_edges = size(mesh.boundary_edges, 1)

    println("=== Mesh Information ===")
    println("Number of nodes:           $num_nodes")
    println("Number of elements:        $num_elements")
    println("Number of boundary edges:  $num_boundary_edges")
    println()
end

"""
Plot the boundary edges of the mesh with optional styling and tag labels.
"""
function plot_domain(
    mesh::Mesh; 
    tag_styles::Union{Dict{Int, Tuple{Symbol, Symbol}}, Nothing} = nothing,
    showtag::Bool = false, tag_text_color::Symbol = :black)

    # Unpack node coordinates
    x = mesh.coordinates[:, 1]
    y = mesh.coordinates[:, 2]
    
    # Initialize empty plot
    plot(; grid=false, aspect_ratio=:equal, legend=false)

    # Loop over all boundary edges and corresponding tags
    for (edge, tag) in zip(eachrow(mesh.boundary_edges), mesh.boundary_tags)
        i1, i2 = edge
        if tag_styles === nothing
             color, linestyle = (:black, :solid)
        else
            color, linestyle = get(tag_styles, tag, (:black, :solid))
        end

        # Get edge endpoints
        x1, y1 = x[i1], y[i1]
        x2, y2 = x[i2], y[i2]

        # Draw edge
        plot!(
            [x1, x2], [y1, y2],
            color=color,
            linestyle=linestyle,
            linewidth=3,
        )

        # Annotate tag at edge midpoint
        if showtag
            xm = (x1 + x2) / 2
            ym = (y1 + y2) / 2
            annotate!(xm, ym, text(string(tag), tag_text_color, 10, :center))
        end
    end

    return current()
end

"""
Plot the full triangular mesh.
"""
function plot_mesh(mesh::Mesh; title::String = "Mesh")
    # Unpack coordinates
    x = mesh.coordinates[:, 1]
    y = mesh.coordinates[:, 2]

    # Transpose element connectivity to 3 × n format
    connectivity = transpose(mesh.elements)

    # Plot mesh
    trimesh(x, y, connectivity; 
        title=title, 
        aspect_ratio=:equal, 
        legend=false
    )
end

"""
Plot the scalar nodal solution on the mesh using fill or contour style.
"""
function plot_solution(
    mesh::Mesh,
    u::AbstractVector;
    style::String = "fill")

    # Extract mesh data
    x = mesh.coordinates[:, 1]
    y = mesh.coordinates[:, 2]
    connectivity = transpose(mesh.elements)  # 3 × n format for TriplotRecipes

    # Plot mesh boundaries
    plot_domain(mesh)

    # Plot solution
    if style == "fill"
        tripcolor!(
            x, y, u, connectivity;
            xlabel="x", ylabel="y", title="Solution")
    elseif style == "contour"
        tricontour!(
            x, y, u, connectivity, 20;
            xlabel="x", ylabel="y", title="Solution")
    else
        error("Unknown plot style: choose \"fill\" or \"contour\"")
    end

    return current()
end

end  # module Mesh2D
