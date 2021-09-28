ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Makie
using CairoMakie
using Plots
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

function nice_divergent_levels(c, clim; nlevels=20)
    levels = range(-clim, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([-cmax], levels, [cmax]))
    return (-clim, clim), levels
end

function nice_nondivergent_levels(c, clim; nlevels=20)
    levels = range(0, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([0], levels, [cmax]))
    return (0, clim), levels
end

function nan_solid(y, z, v, seamount)
    Ny, Nz = size(v)
    y2 = reshape(y, Ny, 1)
    z2 = reshape(z, 1, Nz)
    v[seamount.(0, y2,z2)] .= NaN
    return nothing
end

function show_mask(grid)
    print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

grid = RegularRectilinearGrid(size=(128, 64), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

h0, L              = 0.5, 0.25
seamount(x, y, z)  = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))


V = Field(Center, Center, Center, CPU(), grid_with_seamount)
W = Field(Center, Center, Center, CPU(), grid_with_seamount)


yq = LinRange(-1, 1, 32)
zq = LinRange(-1, 0, 16)


V = [-2*(1 - z/(-2 + h0*exp(-y^2/L^2)))/(-2 + h0*exp(-y^2/L^2)) for y in yq, z in zq]
W = [4*h0*z*y*exp(-y^2/L^2)*(1 - z/(-2 + h0*exp(-y^2/L^2)))/((L^2)*(-2 + h0*exp(-y^2/L^2))^2) for y in yq, z in zq]

mask_immersed_field!(V)
mask_immersed_field!(W)


scene = arrows(yq, zq, V, W, arrowsize = 7, lengthscale = 0.3)
Makie.save("4plot_vel_7_0.3.png", scene)
