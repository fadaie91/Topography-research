ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots
using JLD2
using LinearAlgebra
using Oceananigans.ImmersedBoundaries: solid_node
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

function show_mask(grid)

    print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

grid = RegularRectilinearGrid(size=(16, 8),
                              y=(-1, 1),                        
                              z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

# Gaussian seamount
h0, L = 0.5, 0.25                                        
seamount(x, y, z) = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

#x, y, z, cc = show_mask(grid_with_seamount)
#plt = contourf(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
#savefig(plt,"mask.png")


### Here we define psi and then fill the halo reginos
Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)
h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)

fill_halo_regions!(Ψ, CPU())
### We mask psi
mask_immersed_field!(Ψ)
xs, ys, zs = nodes((Center, Face, Face), grid)
plt = contourf(ys, zs, interior(Ψ)[1,:,:]', xlabel = "y", ylabel = "z", title = "Stream Function")
savefig(plt,"stream-function.png")

#x, y, z, mask = show_mask(grid_with_seamount)
#Ψlims, Ψlevels = nice_nondivergent_levels(Ψ, 1)
#plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Stream Function", color = :jet)
#contour!(plt, ys, zs, interior(Ψ)[1,:,:]', title="Ψ", xlabel="y", ylabel="z",
 #        levels = Ψlevels, clim=Ψlims, color=:black, linewidth=1.5)
#savefig(plt, "Psi.png")

