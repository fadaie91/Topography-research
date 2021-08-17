ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary

grid = RegularRectilinearGrid(size=(256, 128), x=(-12.5e3, 12.5e3), z=(0, 4e3), topology=(Periodic, Flat, Bounded))

# Gaussian bump of width "1"
#H0=0.45
#L=2.5
#H=0.4050
bump(x, y, z) = z < -0.45 + 0.405e3*exp(-(x^2)/(2.5e3^2))

grid_with_bump = ImmersedBoundaryGrid(grid, GridFittedBoundary(bump))

# Tidal forcing
tidal_forcing(x, y, z, t) = 0.0



model = HydrostaticFreeSurfaceModel(architecture = GPU(),
                                    grid = grid_with_bump,
                                    momentum_advection = CenteredSecondOrder(),
                                    free_surface = ImplicitFreeSurface(),
                                    closure = IsotropicDiffusivity(ν=1e-4, κ=1e-5),
                                    tracers = :b,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(f=0),
                                    forcing = (u = tidal_forcing,))

                            

# Linear stratification

set!(model, u = (x, y, z) -> 0.1, b = (x, y, z) -> 1)
#set!(model, u = (x, y, z) -> 0.1, b = (x, y, z) -> -0.1 * (z/0.45))
#set!(model, u = (x, y, z) -> 0.1, b = (x, y, z) -> -0.1 * exp(z/0.1))

progress(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|w|: %.2e",
                             100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                             s.model.clock.time, maximum(abs, model.velocities.w))

gravity_wave_speed = sqrt(model.free_surface.gravitational_acceleration * grid.Lz)
Δt = 0.1 * grid.Δx / gravity_wave_speed
              
simulation = Simulation(model, Δt = Δt, stop_time = 1, progress = progress, iteration_interval = 10)

serialize_grid(file, model) = file["serialized/grid"] = model.grid.grid

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers, (η=model.free_surface.η,)),
                                                      schedule = TimeInterval(0.1),
                                                      prefix = "internal_tide",
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

function nan_solid(x, z, u, bump)
    Nx, Nz = size(u)
    x2 = reshape(x, Nx, 1)
    z2 = reshape(z, 1, Nz)
    u[bump.(x2, 0, z2)] .= NaN
    return nothing
end

function visualize_internal_tide_simulation(prefix)

    filename = prefix * ".jld2"
    file = jldopen(filename)

    grid = file["serialized/grid"]

    # bump(x, y, z) = z < exp(-x^2)
    bump(x, y, z) = z < -0.45 + 0.405e3*exp(-(x^2)/(2.5e3^2))

    xu, yu, zu = nodes((Face, Center, Center), grid)
    xw, yw, zw = nodes((Center, Center, Face), grid)
    xb, yb, zb = nodes((Center, Center, Center), grid)
    xd, yd, zd = nodes((Center, Center, Center), grid)

    b₀ = file["timeseries/b/0"][:, 1, :]

    iterations = parse.(Int, keys(file["timeseries/t"]))    

    anim = @animate for (i, iter) in enumerate(iterations)

        @info "Plotting iteration $iter of $(iterations[end])..."

        u = file["timeseries/u/$iter"][:, 1, :]
        w = file["timeseries/w/$iter"][:, 1, :]
        b = file["timeseries/b/$iter"][:, 1, :]
        η = file["timeseries/η/$iter"][:, 1, 1]
        t = file["timeseries/t/$iter"]

        b′ = b .- b₀
	
	bmax = maximum(abs, b′)
	umax = maximum(abs, u)
	wmax = maximum(abs, w)
    ηmax = maximum(abs, η)


        print("Max b = ", bmax, " Max u = ", umax, " Max w = ", wmax," Max η = ", ηmax, "\n")

        nan_solid(xu, zu, u, bump)
        nan_solid(xw, zw, w, bump)
        nan_solid(xb, zb, b′, bump) 
        nan_solid(xd, zd, η, bump)

        u_title = @sprintf("x velocity, t = %.2f", t)

        u_plot = contourf(xu, zu, u'; title = u_title,                  color = :balance, aspectratio = :equal, linewidth = 0, clim = (-umax, umax))
        w_plot = contourf(xw, zw, w'; title = "z velocity",             color = :balance, aspectratio = :equal, linewidth = 0, clim = (-wmax, wmax))
        b_plot = contourf(xb, zb, b′'; title = "buoyancy perturbation", color = :balance, aspectratio = :equal, linewidth = 0, clim = (-bmax, bmax))
        η_plot = contourf(xd, zd, η'; title = "free surface",           color = :balance, aspectratio = :equal, linewidth = 0, clim = (-ηmax, ηmax))

        plot(u_plot, w_plot, b_plot, η_plot, layout = (4, 1), size = (2400, 2400))
      
    end

    mp4(anim, "seamount.mp4", fps = 16)

    close(file)
end

visualize_internal_tide_simulation("internal_tide")
print("Simulation time = ", prettytime((finish_time - start_time)/1e9), "\n")
