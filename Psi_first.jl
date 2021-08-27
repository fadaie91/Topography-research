ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots

using Oceananigans.ImmersedBoundaries: mask_immersed_field!

function meshgrid(x, y)
    x = reshape(x, 1, length(x))
    y = reshape(y, length(y), 1)
    X = repeat(x, length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

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

x, y, z, cc = show_mask(grid_with_seamount)
plt = contourf(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
savefig(plt,"mask.png")

Psi = Field(Center, Face, Face, CPU(), grid_with_seamount)



Y, Z = meshgrid(y, z)
Psi    = (1 .- Z./(h0*exp.(-Y.^2/L^2) .- 1)).^2



xp, yp, zp = nodes((Center, Face, Face),   grid)

mask_immersed_field!(Psi )

psiplot = contourf(yp, zp, Psi; xlabel = "y", ylabel = "z", title = "Psi")
savefig(psiplot,"psi.png")



