ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots
using JLD2
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

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)

h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
fill_halo_regions!(Ψ, CPU())
mask_immersed_field!(Ψ)

V = YFaceField(CPU(), grid_with_seamount)
W = ZFaceField(CPU(), grid_with_seamount)

set!(V, (x, y, z) -> -2*(1 - z/(-2 + h0*exp(-y^2/L^2)))/(-2 + h0*exp(-y^2/L^2)))
set!(W, (x, y, z) -> 4*h0*z*y*exp(-y^2/L^2)*(1 - z/(-2 + h0*exp(-y^2/L^2)))/((L^2)*(-2 + h0*exp(-y^2/L^2))^2))

fill_halo_regions!(V, CPU())
fill_halo_regions!(W, CPU())
        
mask_immersed_field!(V)
mask_immersed_field!(W)




velocities = PrescribedVelocityFields( v=V, w=W)

model = HydrostaticFreeSurfaceModel(architecture = CPU(),                          
                                    grid = grid_with_seamount,
                                    tracers = :b,
                                   velocities = velocities,
                                    buoyancy = nothing
                                    )



  
B = Field(Center, Center, Center, CPU(), grid_with_seamount)     
set!(B, (x, y, z) -> (1.0 + z) ) 
mask_immersed_field!(B)
xb, yb, zb = nodes((Center, Center, Center), grid)
bplot = contourf(yb, zb, interior(B)[1,:,:]', title="tracer", xlabel="y", ylabel="z")
savefig(bplot, "tracer.png")

Δt = 0.001
set!(model, b = B)                       

progress(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|b|: %.2e",
                             100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                             s.model.clock.time, maximum(abs, model.tracers.b))

xb, yb, zb = nodes((Center, Center, Center), grid)
bplot = contourf(yb, zb, interior(model.tracers.b)[1, :, :]', title="tracer", xlabel="y", ylabel="z")
savefig(bplot, "tracer_initial.png")

simulation = Simulation(model, Δt = Δt, stop_time = 0.3, progress = progress, iteration_interval = 10)

serialize_grid(file, model) = file["serialized/grid"] = model.grid.grid
simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers,
                                                      schedule = TimeInterval(0.02),
                                                      prefix = "flow_over_seamount",
                                                      init = serialize_grid,
                                                      force = true)
                        
start_time = time_ns()
run!(simulation)
finish_time = time_ns()

xb, yb, zb = nodes((Center, Center, Center), grid)
bplot = contourf(yb, zb, interior(model.tracers.b)[1, :, :]'; title="tracer", xlabel="y", ylabel="z", label=@sprintf("t = %.3f", model.clock.time))
savefig(bplot, "tracer_final.png")


@info """
    Simulation complete.
    Output: $(abspath(simulation.output_writers[:fields].filepath))
"""

using JLD2
using Plots
ENV["GKSwstype"] = "100"

function nice_divergent_levels(c, clim; nlevels=20)
    levels = range(-clim, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([-cmax], levels, [cmax]))
    return (-clim, clim), levels
end

function nan_solid(y, z, v, seamount)
    Ny, Nz = size(v)
    y2 = reshape(y, Ny, 1)
    z2 = reshape(z, 1, Nz)
    v[seamount.(0, y2,z2)] .= NaN
    return nothing
end

function visualize_flow_over_seamount_simulation(prefix)

    h0, L = 0.5, 0.25                                        
    seamount(x, y, z) = z < - 1 + h0*exp(-y^2/L^2)
    
    filename = prefix * ".jld2"
    file = jldopen(filename)

    grid = file["serialized/grid"]


    xb, yb, zb = nodes((Center, Center, Center), grid)

    b₀ = file["timeseries/b/0"][1, :, :]

    iterations = parse.(Int, keys(file["timeseries/t"]))    

    anim = @animate for (i, iter) in enumerate(iterations)

        @info "Plotting iteration $iter of $(iterations[end])..."

       
        b = file["timeseries/b/$iter"][1, :, :]
        t = file["timeseries/t/$iter"]

        b′ = b .- b₀
	
	
	bmax = maximum(abs, b′)

        print("Max b = ", bmax)

        
        nan_solid(yb, zb, b′, seamount) 

        v_title = @sprintf("v, t = %.2f", t)

        #b_plot = contourf(yb, zb, b′'; title = "buoyancy perturbation", color = :balance, linewidth = 0, clim = (-bmax, bmax))
	#b_plot = contourf(yb, zb, b′'; title = "buoyancy perturbation")
		b_plot = contourf(yb, zb, b'; title = "buoyancy perturbation")
	# plot(b_plot, layout = (1, 1), size = (1200, 1200))
    end

    mp4(anim, "flow_over_seamount.mp4", fps = 16)

    close(file)
end

visualize_flow_over_seamount_simulation("flow_over_seamount")
print("Simulation time = ", prettytime((finish_time - start_time)/1e9), "\n")
