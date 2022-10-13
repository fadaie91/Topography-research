using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
using Printf
using GLMakie
using Oceananigans.ImmersedBoundaries: mask_immersed_field!


arch = CPU()
#tracer_advection = CenteredSecondOrder()
#momentum_advection = CenteredSecondOrder()
momentum_advection = WENO5()
tracer_advection = WENO5()
#momentum_advection = CenteredFourthOrder()
#tracer_advection = CenteredFourthOrder()
#momentum_advection = UpwindBiasedFirstOrder()
#tracer_advection = UpwindBiasedFirstOrder()
#momentum_advection = UpwindBiasedThirdOrder()
#tracer_advection = UpwindBiasedThirdOrder()
#momentum_advection = UpwindBiasedFifthOrder()
#tracer_advection = UpwindBiasedFifthOrder()


underlying_grid = RectilinearGrid(arch,
                                  size=(16, 8), halo=(3, 3), 
                                  y = (-1, 1),
                                  z = (-1, 0),
                                  topology=(Flat, Periodic, Bounded))

# A bump
h₀ = 0.5 # bump height
L = 0.25 # bump width
@inline h(y) = h₀ * exp(- y^2 / L^2)
@inline seamount(x, y) = - 1 + h(y)

seamount_field = Field{Center, Center, Nothing}(underlying_grid)
set!(seamount_field, seamount)
fill_halo_regions!(seamount_field)

minimum_fractional_Δz = 0.2

immersed_boundaries = [
                       PartialCellBottom(seamount_field.data;
                                         minimum_fractional_Δz),
                       GridFittedBottom(seamount_field.data)
                      ]

b = []
v = []

function progress(sim)
    vmax = maximum(abs, sim.model.velocities.v)
    @info @sprintf("Iter: %d, time: %.2e, max|v|: %.2e",
                   iteration(sim), time(sim), vmax)

    return nothing
end

tracer_errors  = Dict()

#for ib in immersed_boundaries
    #grid = ImmersedBoundaryGrid(underlying_grid, ib)
    #grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data))
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(seamount_field.data))

    @show grid
    ### In this example of the paper the velocities were prescribe
    ### which means that they won't evolve during the simulation
    ### and using the function below we prescribe the velocities
    B= Field{Center, Center, Center}(grid)   
    N² = 1
    ###I should go back to adcroft==> see what they prescribed for bouyancy                                  
    set!(B, (x, y, z) -> (N² * z) ) 
    mask_immersed_field!(B)

    #U, V, W = model.velocities

    Ψ = Field{Center, Face, Face}(grid)

    h(y)    = h₀*exp(-y^2/L^2)
    ζ(y, z) = z/(h(y) - 1)
    set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
    fill_halo_regions!(Ψ, arch)


    ### We mask psi
    mask_immersed_field!(Ψ)

    ### V is (Center, Face, Center)
    V = YFaceField(grid)
    ### W is (Center, Centere, Face)
    W = ZFaceField(grid)

    ### V and W are deravatives of psi 
    ### which in this code we use exact expression 
    ### which we computed in maple
    ### we are calling this method 'analytical'
    V.=  ∂z(Ψ)
    W.= -∂y(Ψ)

    ### We mask V and W        
    mask_immersed_field!(V)
    mask_immersed_field!(W)

    ### We fill the halo regions of V and W

    fill_halo_regions!(V, arch)
    fill_halo_regions!(W, arch)

    velocities = PrescribedVelocityFields(v=V, w=W)
    
    model = HydrostaticFreeSurfaceModel(; grid,
                                            tracer_advection,
                                            momentum_advection,
                                            #coriolis = FPlane(f=0.1),
                                            tracers = :b,
                                            velocities = velocities,
                                            buoyancy = BuoyancyTracer())


    
    #set!(model, b = B,  w=W, v=V)
    set!(model, b = B)
    tracer_initial= sum(interior(model.tracers.b))

    simulation = Simulation(model; Δt=1e-3, stop_time=1)
    simulation.callbacks[:p] = Callback(progress, IterationInterval(10))


    run!(simulation)
