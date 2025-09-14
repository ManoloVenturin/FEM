# basic_mesh.jl

coordinates = [
    0.0  0.0;
    1.0  0.0;
    0.0  1.0;
    1.0  1.0
]

elements = [
    1  2  3;
    2  4  3
]

boundaries = [
    1 2 "Dirichlet";
    2 4 "Neumann";
    4 3 "Neumann";
    3 1 "Dirichlet"
]


# %%


using LinearAlgebra
using Plots, TriplotRecipes

# Unpack data
x = coordinates[:,1]
y = coordinates[:,2]

# Mesh connectivity require a format 3 x n_elem
t = transpose(elements)

# Show mesh
trimesh(x, y, t; title="Mesh", aspect_ratio=:equal, legend=false)

# Save fig
filename = joinpath(@__DIR__, "img", "basic_mesh.png")
savefig(filename)


# %%


"""
Plot domain from boundaries
"""
function plot_domain(x, y, boundaries; mapbcs=true)
    style_map = Dict(
        "Dirichlet" => (:red, :solid),
        "Neumann" => (:blue, :dash)
    )

    plot(; grid=false, aspect_ratio=:equal)
    
    for row in eachrow(boundaries)
        i1 = row[1]
        i2 = row[2]
        cond = row[3]
        color, linestyle =  (:black, :solid)
        if mapbcs
            color, linestyle = get(style_map, cond, (:black, :solid))
        end
        plot!([x[i1], x[i2]], [y[i1], y[i2]], color=color, 
            linestyle=linestyle, linewidth=4, label=false)
    end
    return current()
end

# Show boundary and save it
plot_domain(x, y, boundaries)
filename = joinpath(@__DIR__, "img", "basic_boundary.png")
savefig(filename)


# %%


# Nodal solution
z = [0.; 1.; 2; 3]

# Solution visualization (fill/heatmap)
plot_domain(x, y, boundaries; mapbcs=false)
tripcolor!(x, y, z, t; xlabel="x", ylabel="y", title="Solution")

# Save fig
filename = joinpath(@__DIR__, "img", "basic_plot1.png")
savefig(filename)


# Solution visualization (contour lines)
plot_domain(x, y, boundaries; mapbcs=false)
tricontour!(x, y, z, t, 20; xlabel="x", ylabel="y", title="Solution")

# Save fig
filename = joinpath(@__DIR__, "img", "basic_plot2.png")
savefig(filename)
