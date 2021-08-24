#ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary

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

U(x, y, z) =  0.0
V(x, y, z) = -2*(1-z/(-2+h0*exp(-y^2/L^2)))/(-2+h0*exp(-y^2/L^2))
W(x, y, z) =  4*(1-z/(-2+h0*exp(-y^2/L^2)))/(-2+h0*exp(-y^2/L^2))^2*y*z*h0/L^2*exp(-y^2/L^2)

velocities = PrescribedVelocityFields(u=U, v=V, w=W)

model = HydrostaticFreeSurfaceModel(architecture = CPU(),                          
                                    grid = grid_with_seamount,
                                    momentum_advection = CenteredSecondOrder(), 
                                    free_surface = ImplicitFreeSurface(),
                                    closure = IsotropicDiffusivity(ν=1e-4, κ=0),
                                    tracers = :b,
#                                    velocities = velocities,
                                    buoyancy = BuoyancyTracer())

g = model.free_surface.gravitational_acceleration
gravity_wave_speed = sqrt(g * grid.Lz)
Δt = 0.01 * grid.Δy / gravity_wave_speed              

set!(model, u = U, v = V, w = W, b = 1.0)                       
progress(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|w|: %.2e",
                             100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                             s.model.clock.time, maximum(abs, model.velocities.w))

simulation = Simulation(model, Δt = Δt, stop_time = 0.2, progress = progress, iteration_interval = 10)

serialize_grid(file, model) = file["serialized/grid"] = model.grid.grid

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers, (η=model.free_surface.η,)),
                                                      schedule = TimeInterval(0.02),
                                                      prefix = "flow_over_seamount",
                                                      init = serialize_grid,
                                                      force = true)
                        
start_time = time_ns()
run!(simulation)
finish_time = time_ns()

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

    xv, yv, zv = nodes((Face, Center, Center), grid)
    xw, yw, zw = nodes((Center, Center, Face), grid)
    xb, yb, zb = nodes((Center, Center, Center), grid)

    b₀ = file["timeseries/b/0"][1, :, :]

    iterations = parse.(Int, keys(file["timeseries/t"]))    

    anim = @animate for (i, iter) in enumerate(iterations)

        @info "Plotting iteration $iter of $(iterations[end])..."

        v = file["timeseries/v/$iter"][1, :, :]
        w = file["timeseries/w/$iter"][1, :, :]
        η = file["timeseries/η/$iter"][1, :, 1]
        b = file["timeseries/b/$iter"][1, :, :]
        t = file["timeseries/t/$iter"]

        b′ = b .- b₀
	
	vmax = maximum(abs, v)
	wmax = maximum(abs, w)
	bmax = maximum(abs, b′)

        print("Max b = ", bmax, " Max v = ", vmax, " Max w = ", wmax, "\n")

        nan_solid(yv, zv, v, seamount)
        nan_solid(yw, zw, w, seamount)
        nan_solid(yb, zb, b′, seamount) 

        v_title = @sprintf("v, t = %.2f", t)

        v_plot = contourf(yv, zv, v'; title = v_title,                 color = :balance, linewidth = 0, clim = (-vmax, vmax))
        w_plot = contourf(yw, zw, w'; title = "w",                     color = :balance, linewidth = 0, clim = (-wmax, wmax))
        b_plot = contourf(yb, zb, b′'; title = "buoyancy perturbation", color = :balance, linewidth = 0, clim = (-bmax, bmax))

        plot(v_plot, w_plot, b_plot, layout = (3, 1), size = (1200, 1200))
    end

    mp4(anim, "flow_over_seamount.mp4", fps = 16)

    close(file)
end

visualize_flow_over_seamount_simulation("flow_over_seamount")
print("Simulation time = ", prettytime((finish_time - start_time)/1e9), "\n")

