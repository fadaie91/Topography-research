using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.ImmersedBoundaries: solid_node

grid = RegularRectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

h0, L              = 0.5, 0.25
seamount(x, y, z)  = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)

h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)

mask_immersed_field!(Ψ)
#rotl90=>transpose the psi
rotl90(interior(Ψ)[1,:,:])

#when solid_node is true it means that it is inside the seamount
solid_node(Center(), Face(), Face(), 1,6,2, grid_with_seamount)
#when solid_node is true it means that it is outside the seamount
solid_node(Center(), Face(), Face(), 1,6,3, grid_with_seamount)
