ENV["GKSwstype"] = "nul"
using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.BoundaryConditions: fill_halo_regions!

using Plots

using Oceananigans.ImmersedBoundaries: mask_immersed_field!

function show_mask(grid)

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
plt = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
savefig(plt,"mask.png")

h(y) = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
fill_halo_regions!(Ψ, CPU())
V = YFaceField(CPU(), grid_with_seamount)
W = ZFaceField(CPU(), grid_with_seamount)
V.=  ∂z(Ψ)
W.= -∂y(Ψ)

xv, yv, zv = nodes((Face, Center, Center), grid)
xw, yw, zw = nodes((Center, Center, Face), grid)
xp, yp, zp = nodes((Center, Face, Face),   grid)

mask_immersed_field!(Ψ)
mask_immersed_field!(V)
mask_immersed_field!(W)

plt = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region", color = :jet)
contour!(plt, yp, zp, interior(Ψ)[1,:,:]', title="Ψ", xlabel="y", ylabel="z", clim=(0, 1), color=:black, linewidth=1.5)
savefig(plt, "Psi.png")

plt = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region", color = :jet)
contour!(plt, yv, zv, interior(V)[1,:,:]', title="V", xlabel="y", ylabel="z", color=:black, linewidth=1.5)
savefig(plt, "V.png")

plt = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region", color = :jet)
contour!(plt, yw, zw, interior(W)[1,:,:]', title="W", xlabel="y", ylabel="z", color=:black, linewidth=1.5)
savefig(plt, "W.png")

velocities = PrescribedVelocityFields(v=V, w=W)

model = HydrostaticFreeSurfaceModel(architecture = CPU(),                          
                                    grid = grid_with_seamount,
                                    momentum_advection = CenteredSecondOrder(), 
                                    #free_surface = ImplicitFreeSurface(),
                                    closure = IsotropicDiffusivity(ν=1e-4, κ=0),
                                    tracers = :b,
                                    velocities = velocities,
                                    buoyancy = BuoyancyTracer())


Δt = 0.001

set!(model, b = (x, y, z) -> 1 + z)                       
progress(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|w|: %.2e",
                             100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                             s.model.clock.time, maximum(abs, model.velocities.w))

simulation = Simulation(model, Δt = Δt, stop_time = 0.3, progress = progress, iteration_interval = 10)

pop!(simulation.diagnostics)

serialize_grid(file, model) = file["serialized/grid"] = model.grid.grid

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers,
                                                      schedule = TimeInterval(0.005),
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

        b = file["timeseries/b/$iter"][1, :, :]
        t = file["timeseries/t/$iter"]

        b′ = b .- b₀
	
        bmax = maximum(abs, b′)

        print("Max b = ", bmax, "\n")

        nan_solid(yb, zb, b′, seamount) 

        b_plot = contour(yb, zb, b′'; title = "buoyancy perturbation", linewidth = 1, clim = (-bmax, bmax))
    end

    mp4(anim, "flow_over_seamount.mp4", fps = 16)

    close(file)
end

visualize_flow_over_seamount_simulation("flow_over_seamount")
print("Simulation time = ", prettytime((finish_time - start_time)/1e9), "\n")

