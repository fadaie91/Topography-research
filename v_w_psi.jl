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

grid = RegularRectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

h0, L              = 0, 0.25
seamount(x, y, z)  = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

x, y, z, mask = show_mask(grid_with_seamount)
plt = heatmap(y, z, interior(mask)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region")
savefig(plt,"mask.png")

Ψ = Field(Center, Face, Face, CPU(), grid_with_seamount)
V = YFaceField(CPU(), grid_with_seamount)
#W = ZFaceField(CPU(), grid_with_seamount)
W = Field(Center, Face, Center, CPU(), grid_with_seamount)
fill_halo_regions!(Ψ, CPU())
fill_halo_regions!(V, CPU())
fill_halo_regions!(W, CPU())
h(y)    = h0*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
#fill_halo_regions!(Ψ, CPU())

#(Ψ)[1,:,:]'
#interior(Ψ)[1,:,:]'


set!(V, (x, y, z) -> -2*(1 - z/(-2 + h0*exp(-y^2/L^2)))/(-2 + h0*exp(-y^2/L^2)))

set!(W, (x, y, z) -> 4*h0*z*y*exp(-y^2/L^2)*(1 - z/(-2 + h0*exp(-y^2/L^2)))/((L^2)*(-2 + h0*exp(-y^2/L^2))^2))

#set!(W, (x, y, z) -> 4*h0*z*y*exp(-y^2/L^2))

#fill_halo_regions!(V, CPU())
#fill_halo_regions!(W, CPU())

xp, yp, zp = nodes((Center, Face, Face),   grid)
xv, yv, zv = nodes((Face, Center, Center), grid)
#xv, yv, zv = nodes((Center, Face, Face), grid)
xw, yw, zw = nodes((Center, Face, Center), grid)
#xw, yw, zw = nodes((Center, Face, Face), grid)


mask_immersed_field!(Ψ)
mask_immersed_field!(V)
mask_immersed_field!(W)

psiplot = contourf(yp, zp, interior(Ψ)[1,:,:]', title="Psi_numerical", xlabel="y", ylabel="z")
savefig(psiplot, "psinumerical.png")


plt_v = quiver(yv, zv, interior(V)[1,:,:]', title="V_numerical")
savefig(plt_v, "quivervnumeical.png")

plt_w = quiver(yw, zw, interior(W)[1,:,:]', title="w_numerical")
savefig(plt_w, "quiverwnumeical.png")


(W)[1,:,:]'
interior(W)[1,:,:]'


