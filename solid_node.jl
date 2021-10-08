using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.ImmersedBoundaries: solid_node
using Oceananigans.BoundaryConditions: fill_halo_regions!

grid = RegularRectilinearGrid(size=(8, 4), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

h0, L              = 0.5, 0.25
seamount(x, y, z)  = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)

h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
fill_halo_regions!(Ψ, CPU())
  
mask_immersed_field!(Ψ)
#mask_immersed_field!(interior(Ψ), NaN)
#rotl90=>transpose the psi
rotl90(interior(Ψ)[1,:,:])


V = YFaceField(CPU(), grid_with_seamount)
W = ZFaceField(CPU(), grid_with_seamount)

V.=  ∂z(Ψ)
W.= -∂y(Ψ)

fill_halo_regions!(V, CPU())
fill_halo_regions!(W, CPU())

interior(V)[1,:,:]
interior(W)[1,:,:]

D = Field(Center, Center, Center, CPU(), grid_with_seamount)
D .= ∂y(V) + ∂z(W)

D1 = Field(Center, Center, Center, CPU(), grid_with_seamount)
D1 .= ∂y(V) 

D2 = Field(Center, Center, Center, CPU(), grid_with_seamount)
D2 .=  ∂z(W)
interior(D)[1,:,:]
#when solid_node is true it means that it is inside the seamount
solid_node(Center(), Face(), Face(), 1,6,2, grid_with_seamount)
#when solid_node is true it means that it is outside the seamount
solid_node(Center(), Face(), Face(), 1,6,3, grid_with_seamount)
