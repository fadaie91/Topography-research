#ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots
using JLD2
using LinearAlgebra
using Oceananigans.ImmersedBoundaries: solid_node
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

include("needed_functions.jl")

function nice_nondivergent_levels(c, clim; nlevels=20)
    levels = range(0, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([0], levels, [cmax]))
    return (0, clim), levels
end



### Define grid with a seamount
h0, L = 0.5, 0.25
seamount(x, y, z) = z < - 1 + h0*exp(-y^2/L^2)
grid = RectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),
                       topology=(Flat, Periodic, Bounded))
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))
x, y, z, mask = show_mask(grid_with_seamount)
### Plot masked region: mask.png
x, y, z, cc = show_mask(grid_with_seamount)
plt = contourf(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z",
                title = "Masked Region")
savefig(plt,"mask.png")

### Define Fields
Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)
V = YFaceField(CPU(), grid_with_seamount)
W = ZFaceField(CPU(), grid_with_seamount)
B = Field(Center, Center, Center, CPU(), grid_with_seamount)

h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)

set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
#V.=  ∂z(Ψ)
#W.= -∂y(Ψ)
set!(V, (x, y, z) -> -2 * (1 - ( z / (-1+h0*exp(-y^2/L^2)) ) )/(-1+h0*exp(-y^2/L^2)))
set!(W, (x, y, z) -> 4*h0*z*y*exp(-y^2/L^2)*(1 - z/(-1 + h0*exp(-y^2/L^2)))/((L^2)*(-1 + h0*exp(-y^2/L^2))^2))
set!(B, (x, y, z) -> (1.0 + z) )

fill_halo_regions!(Ψ, CPU())
fill_halo_regions!(V, CPU())
fill_halo_regions!(W, CPU())
fill_halo_regions!(B, CPU()) 

mask_immersed_field!(Ψ)
mask_immersed_field!(V)
mask_immersed_field!(W)
mask_immersed_field!(B)

Ψlims, Ψlevels = nice_nondivergent_levels(Ψ, 1)
xp, yp, zp = nodes((Center, Face, Face),   grid)

plt =heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Psi", color = :jet)
contour!(plt, yp, zp, interior(Ψ)[1,:,:]', title="Ψ", xlabel="y", ylabel="z") 


savefig(plt, "Psi.png")
