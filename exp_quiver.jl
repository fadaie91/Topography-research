ENV["GKSwstype"] = "nul"
using Plots

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
        
vectorfield = vectorfield2d
function meshgrid(n)
    xs = ones(n) .* (1:n)'
    ys = xs'
    xys = permutedims(cat(xs, ys; dims = 3), [3, 1, 2])
    return reshape(xys, 2, n^2)
end
h0, L              = 0, 0.25
vf(y, z) = -2*(1 - z/(-2 + h0*exp(-y^2/L^2)))/(-2 + h0*exp(-y^2/L^2)) # circular vectorfield

grid = RegularRectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))
exp_quiver=vectorfield2d(vf, grid)
savefig(exp_quiver, "exp_quiver.png")