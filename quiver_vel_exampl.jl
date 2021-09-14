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
function vectorfield2d(field, points, arrowlength=0.1)
    # More input pattern parsing is solved by the Plots package, but I don't know how.
    errormessage = "Incorrect formatting of points. Please format them as [x1 y1; x2, y2;...]"
    
    if typeof(points) <: Array{<:Number, 2} && size(points)[1] === 2
        vectors = similar(points)
        for i in 1:size(points)[2]
            vectors[:, i] .= collect(field(points[:, i]...))
        end
    else
        error(errormessage)
    end
    vectors .*= arrowlength
    quiver(points[1, :],points[2, :],quiver=(vectors[1, :], vectors[2, :]))
   # display(plot!())
end

grid = RegularRectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

h0, L              = 0.1, 0.25
seamount(x, y, z)  = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

x, y, z, mask = show_mask(grid_with_seamount)
plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
savefig(plt,"mask.png")

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)
V = YFaceField(CPU(), grid_with_seamount)
W = Field(Center, Face, Center, CPU(), grid_with_seamount)

fill_halo_regions!(Ψ, CPU())
fill_halo_regions!(V, CPU())
fill_halo_regions!(W, CPU())
            
h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)



meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

y, z = meshgrid(-1:1/16:1, -1:1/8:0)

V = @. -2*(1 - z/(-2 + h0*exp(-y^2/L^2)))/(-2 + h0*exp(-y^2/L^2))
W = @. 4*h0*z*y*exp(-y^2/L^2)*(1 - z/(-2 + h0*exp(-y^2/L^2)))/((L^2)*(-2 + h0*exp(-y^2/L^2))^2)




xt, yt, zt = nodes((Center, Center, Center), grid)


plt_vel=quiver(yt, zt, quiver=(V, W))
savefig(plt_vel, "quivervelnumeical.png")


plt_v=quiver(yt, zt, quiver=V)
savefig(plt_v, "plt_v.png")

plt_w=quiver(yt, zt, quiver=W)
savefig(plt_w, "plt_w.png")
