ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary

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

Y, Z = meshgrid(y, z)
Psi    = (1 .- Z./(h0*exp.(-Y.^2/L^2) .- 1)).^2

psiplot = contourf(y, z, Psi; xlabel = "y", ylabel = "z", title = "Psi")
savefig(psiplot,"psi.png")



Y, Z = meshgrid(y, z)
V    = -2*(1 .- Z./(-2 .+ h0*exp.(-Y.^2/L^2)))./(-2 .+ h0*exp.(-Y.^2/L^2))



vplot = contourf(y, z, V; xlabel = "y", ylabel = "z", title = "V")
savefig(vplot,"v.png")



Y, Z = meshgrid(y, z)
W    = 4*h0.*Z.*Y.*exp.(-Y.^2/L^2).*(1 .- Z./(-2 .+ h0*exp.(-Y.^2/L^2)))./((L^2)*(-2 .+ h0*exp.(-Y.^2/L^2)).^2)
#W(x, y, z) =  4*(1-z/(-2+h0*exp(-y^2/L^2)))/(-2+h0*exp(-y^2/L^2))^2*y*z*h0/L^2*exp(-y^2/L^2)

wplot = contourf(y, z, W; xlabel = "y", ylabel = "z", title = "W")
savefig(wplot,"w.png")