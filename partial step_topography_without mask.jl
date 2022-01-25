using Oceananigans
using Plots


h0, L = 0.5, 0.25
topography(y) = -1 + h0*exp(-y^2/L^2)
grid = RectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),
                                     topology=(Flat, Periodic, Bounded), halo=(3,3))

height = (topography.(grid.yᵃᶠᵃ[1:grid.Ny]) + 4*topography.(grid.yᵃᶜᵃ[1:grid.Ny]) + topography.(grid.yᵃᶠᵃ[2:grid.Ny+1]))/6        
plt=plot(grid.yᵃᶜᵃ[1:grid.Ny], height, lw=4)
