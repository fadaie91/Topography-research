ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!

using Plots

using Oceananigans.ImmersedBoundaries: mask_immersed_field!

function nice_divergent_levels(c, clim; nlevels=20)
    levels = range(-clim, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([-cmax], levels, [cmax]))
    return (-clim, clim), levels
end

function nice_nondivergent_levels(c, clim; nlevels=20)
    levels = range(0, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([0], levels, [cmax]))
    return (0, clim), levels
end

function nan_solid(y, z, v, seamount)
    Ny, Nz = size(v)
    y2 = reshape(y, Ny, 1)
    z2 = reshape(z, 1, Nz)
    v[seamount.(0, y2,z2)] .= NaN
    return nothing
end

function show_mask(grid)
    print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

grid = RegularRectilinearGrid(size=(128, 64), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

h0, L              = 0.5, 0.25
seamount(x, y, z)  = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

x, y, z, mask = show_mask(grid_with_seamount)
plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
savefig(plt,"mask.png")

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)
V = YFaceField(CPU(), grid_with_seamount)
W = Field(Center, Center, Face, CPU(), grid_with_seamount)

h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
set!(V, (x, y, z) -> -2*(1 - z/(-2 + h0*exp(-y^2/L^2)))/(-2 + h0*exp(-y^2/L^2)))
set!(W, (x, y, z) -> 4*h0*z*y*exp(-y^2/L^2)*(1 - z/(-2 + h0*exp(-y^2/L^2)))/((L^2)*(-2 + h0*exp(-y^2/L^2))^2))



xp, yp, zp = nodes((Center, Face, Face),   grid)
xv, yv, zv = nodes((Center, Face, Center), grid)
xw, yw, zw = nodes((Center, Center, Face), grid)


mask_immersed_field!(Ψ)
mask_immersed_field!(V)
mask_immersed_field!(W)

Ψlims, Ψlevels = nice_nondivergent_levels(Ψ, 1)

plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region", color = :jet)
contour!(plt, yp, zp, interior(Ψ)[1,:,:]', title="Ψ", xlabel="y", ylabel="z",
         levels = Ψlevels, clim=Ψlims, color=:black, linewidth=1.5)
savefig(plt, "Psi.png")

plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region", color = :jet)
contour!(plt, yv, zv, interior(V)[1,:,:]', title="V", xlabel="y", ylabel="z", color=:black, linewidth=1.5)
savefig(plt, "V.png")

plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region", color = :jet)
contour!(plt, yw, zw, interior(W)[1,:,:]', title="W", xlabel="y", ylabel="z", color=:black, linewidth=1.5)
savefig(plt, "W.png")

velocities = PrescribedVelocityFields(v=V, w=W)

model = HydrostaticFreeSurfaceModel(architecture = CPU(),                          
                                    grid = grid_with_seamount,
                                    tracers = :θ,
                                    velocities = velocities,
                                    buoyancy = nothing
                                    )

Δt = 0.001

set!(model, θ = (x, y, z) -> 1 + z)                       
progress(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|b|: %.2e",
                             100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                             s.model.clock.time, maximum(abs, model.tracers.θ))

simulation = Simulation(model, Δt = Δt, stop_time = 1, progress = progress, iteration_interval = 10)

pop!(simulation.diagnostics)

serialize_grid(file, model) = file["serialized/grid"] = model.grid.grid

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers,
                                                      schedule = TimeInterval(0.01),
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

function visualize_flow_over_seamount_simulation(prefix)

    h0, L = 0.5, 0.25                                        
    seamount(x, y, z) = z < - 1 + h0*exp(-y^2/L^2)
    
    filename = prefix * ".jld2"
    file = jldopen(filename)

    grid = file["serialized/grid"]

    xv, yv, zv = nodes((Center, Face, Center), grid)
    xw, yw, zw = nodes((Center, Center, Face), grid)
    xθ, yθ, zθ = nodes((Center, Center, Center), grid)

    iterations = parse.(Int, keys(file["timeseries/t"]))    

    anim = @animate for (i, iter) in enumerate(iterations)

        @info "Plotting iteration $iter of $(iterations[end])..."

        θ = file["timeseries/θ/$iter"][1, :, :]
        t = file["timeseries/t/$iter"]

        θmax = maximum(θ)
        θmin = minimum(θ)
        print("Max θ = ", θmax, " and Min θ = ", θmin, "\n")

        nan_solid(yθ, zθ, θ, seamount) 

        plt = contour(yθ, zθ, θ'; title = "tracer", linewidth = 1)
    end

    mp4(anim, "flow_over_seamount.mp4", fps = 16)

    close(file)
end

visualize_flow_over_seamount_simulation("flow_over_seamount")
print("Simulation time = ", prettytime((finish_time - start_time)/1e9), "\n")
