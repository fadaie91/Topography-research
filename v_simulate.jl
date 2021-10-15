ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots

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

x, y, z, cc = show_mask(grid_with_seamount)
plt = contourf(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
savefig(plt,"mask.png")

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)

h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
fill_halo_regions!(Ψ, CPU())
mask_immersed_field!(Ψ)

V = YFaceField(CPU(), grid_with_seamount)
W = ZFaceField(CPU(), grid_with_seamount)

V.=  ∂z(Ψ)
W.= -∂y(Ψ)

fill_halo_regions!(V, CPU())
fill_halo_regions!(W, CPU())
        
mask_immersed_field!(V)
mask_immersed_field!(W)

xv, yv, zv = nodes((Center, Face, Center), grid)
vplot = contourf(yv, zv, interior(V)[1,:,:]', title="v", xlabel="y", ylabel="z")
savefig(vplot, "v.png")


velocities = PrescribedVelocityFields( v=V, w=W)

model = HydrostaticFreeSurfaceModel(architecture = CPU(),                          
                                    grid = grid_with_seamount,
                                    tracers = :θ,
                           #         velocities = velocities, is there still a problem with prescribed velocity???
                                    buoyancy = nothing
                                    )



Δt = 0.001             

set!(model,  v = V, w = W)                       
progress(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|v|: %.2e",
                             100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                             s.model.clock.time, maximum(abs, model.velocities.v))

simulation = Simulation(model, Δt = Δt, stop_time = 0.2, progress = progress, iteration_interval = 10)

serialize_grid(file, model) = file["serialized/grid"] = model.grid.grid

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                                                      schedule = TimeInterval(0.02),
                                                      prefix = "flow_over_seamount",
                                                      init = serialize_grid,
                                                      force = true)
                        
start_time = time_ns()
run!(simulation)
finish_time = time_ns()


xv, yv, zv = nodes((Center, Face, Center), grid)
plt_v_final=plot!(vplot, yv ,zv, interior(V)[1,:,:]', linewidth=2,
      label=@sprintf("t = %.3f", model.clock.time))
savefig(plt_v_final, "v_final.png")
