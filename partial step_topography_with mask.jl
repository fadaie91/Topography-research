using Oceananigans
using Plots
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

   
function show_mask(grid)

    #print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

h0, L = 0.5, 0.25
topography(y) = -1 + h0*exp(-y^2/L^2)
grid = RectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),
                                     topology=(Flat, Periodic, Bounded), halo=(3,3))

height = (topography.(grid.yᵃᶠᵃ[1:grid.Ny]) + 4*topography.(grid.yᵃᶜᵃ[1:grid.Ny]) + topography.(grid.yᵃᶠᵃ[2:grid.Ny+1]))/6        
plt=plot(grid.yᵃᶜᵃ[1:grid.Ny], height, lw=4)
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(height))

### Plot masked region: mask.png
x, y, z, c = show_mask(grid_with_seamount)
plt = heatmap(y, z, interior(c)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
