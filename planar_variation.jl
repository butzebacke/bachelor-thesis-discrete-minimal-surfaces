using GLMakie
using Observables
using LinearAlgebra
using Colors
using GeometryBasics

n, m = 5,5
x0, y0 = -2, -2
dist = 1
show_isp = false
fig = Figure()
isaxlabelvisible = false
ax = Axis(fig[1,1], aspect=DataAspect())
ax3d = Axis3(fig[1,2], xlabelvisible = isaxlabelvisible,
    ylabelvisible = isaxlabelvisible,
    zlabelvisible = isaxlabelvisible,
    xticksvisible = isaxlabelvisible,
    yticksvisible = isaxlabelvisible,
    zticksvisible = isaxlabelvisible,
    xticklabelsvisible = isaxlabelvisible,
    yticklabelsvisible = isaxlabelvisible,
    zticklabelsvisible = isaxlabelvisible
   )
ax3d.aspect = :data
display(fig)

Makie.deregister_interaction!(ax, :rectanglezoom) 

function invStereographicProj(p::Point2f)
    x,y = p
    sqnorm = x^2 + y^2
    return Point3f(2*x / (sqnorm + 1), 2*y / (sqnorm + 1), (sqnorm - 1)/(sqnorm + 1))
end

function circle_or_line_coords_from_3pts(p1,p2,p3)
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    # Check collinearity via determinant (area of triangle = 0)
    d = 2 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
    tol = 0.05

    if abs(d) < tol # modify to have less points
        v = p1 - p2
        coords = [p1 - 10*v, p2 + 10*v]
    else
        n = 50

        B = (x1^2 + y1^2)*(y3 - y2) + (x2^2 + y2^2)*(y1 - y3) + (x3^2 + y3^2)*(y2 - y1)
        C = (x1^2 + y1^2)*(x2 - x3) + (x2^2 + y2^2)*(x3 - x1) + (x3^2 + y3^2)*(x1- x2)

        xc = -B / d
        yc = -C / d

        center = xc + im * yc
        radius = norm(center .- (x1 + im * y1))

        coords = Point2f.(reim.([radius * exp(im * 2*pi/n * k) + center for k in 0:n]))
    end
    return coords
end

# calculate p12 s.t. cr(p,p1,p12,p2) = mu, watch out for division by 0
function calculate_mu_p12_pos(p::Point2f, p1::Point2f, p2::Point2f, mu)
    p_complex = Complex(p[1],p[2])
    p1_complex = Complex(p1[1],p1[2])
    p2_complex = Complex(p2[1],p2[2])
    p12_complex = (-(p_complex*p2_complex) + p1_complex*p2_complex + p_complex*p1_complex*mu - p1_complex*p2_complex*mu)/
    (-p_complex + p1_complex + p_complex*mu - p2_complex*mu)
    return Point2f(real(p12_complex),imag(p12_complex))
end

# takes in edges Observable and points matrix and updates edges
function update_edges!(edges, points::AbstractMatrix{P}) where {P<:Point}
    segs = Vector{P}()
    for i in 1:n, j in 1:m
        if i < n
            push!(segs, points[i,j])
            push!(segs, points[i+1,j])
        end
        if j < m
            push!(segs, points[i,j])
            push!(segs, points[i,j+1])
        end
    end
    edges[] = segs
end

# f should be n x m array with Point3f entries
function discrete_dual(f)
    (n, m) = size(f)
    f_star0 = Point3f(0)
    # index [i,j] becomes n*(j-1) + i
    # del_f1[k] = f_i+1,j - f_i,j
    del1_f = Array{Point3f}(undef,n-1,m)
    for i in 1:(n-1), j in 1:m
        del1_f[i,j] = f[i+1,j] - f[i,j] 
    end
    del2_f = Array{Point3f}(undef,n,m-1)
    for i in 1:n, j in 1:(m-1)
        del2_f[i,j] = f[i,j+1] - f[i,j] 
    end
    # second deltas already flipped sign
    del1_f_norm = norm.(del1_f)
    del2_f_norm = norm.(del2_f)
    
    del1_f_star = similar(del1_f)
    del2_f_star = similar(del2_f)

    for i in 1:(n-1), j in 1:m
        del1_f_norm[i,j]==0 ? val = Point3f(0) : val = del1_f[i,j] / (del1_f_norm[i,j]^2)
        del1_f_star[i,j] = val
    end

    for i in 1:n, j in 1:(m-1)
        del2_f_norm==0 ? val = Point3f(0) : val = - del2_f[i,j] / (del2_f_norm[i,j]^2)
        del2_f_star[i,j] = val
    end

    f_star = similar(f)
    f_star[1,1] = f_star0
    
    for j in 2:m
        f_star[1,j] = f_star[1,j-1] + del2_f_star[1,j-1]
    end
    for i in 2:n, j in 1:m
        f_star[i,j] = f_star[i-1,j] + del1_f_star[i-1,j]
    end

    # normalize
    p_avg = reduce(+,f_star) / (n*m)
    f_star = f_star .- p_avg

    return f_star
end

points_matrix = [Point2f(x0 + (i-1)*dist, y0 + (j-1)*dist) for i in 1:n, j in 1:m]
points_matrix_isp = invStereographicProj.(points_matrix)
points_matrix_isp_dual = discrete_dual(points_matrix_isp)

obs_points_matrix = Observable(points_matrix)
obs_points_matrix_isp = Observable(points_matrix_isp)
obs_points_matrix_isp_dual = Observable(points_matrix_isp_dual)

obs_edges = Observable(Vector{Point2f}())
update_edges!(obs_edges, points_matrix)

obs_edges_isp = Observable(Vector{Point3f}())
update_edges!(obs_edges_isp, points_matrix_isp)

obs_edges_isp_dual = Observable(Vector{Point3f}())
update_edges!(obs_edges_isp_dual, points_matrix_isp_dual)

p_move = Observable(Point2f(0))
p_move[] = points_matrix[1,1]

linesegments!(ax, obs_edges, color=:steelblue)

p_move_scat = scatter!(ax, p_move, color="green", markersize=15)

scatter!(ax, lift(x -> vec(x)[3:2:end], obs_points_matrix); color=:white, strokecolor=:black, strokewidth=1, markersize = 6)
scatter!(ax, lift(x -> vec(x)[2:2:end], obs_points_matrix); color=:black, strokecolor=:black, strokewidth=1, markersize = 6)

if show_isp scatter!(ax3d, lift(x -> vec(x), obs_points_matrix_isp)) end
scatter!(ax3d, lift(x -> vec(x), obs_points_matrix_isp_dual))

if show_isp linesegments!(ax3d, obs_edges_isp) end
linesegments!(ax3d, obs_edges_isp_dual)

# Calculate all other white vertices to have cr = mu requirement, combinatorics are kind of complicated here, we kind of calculate too much here
function update_grid_cr_mu!(grid, mu)
    for j in 1:(m-1), i in 1:(n-1)
        if j%2 == 1
            if i%2 == 1
                #i odd, j odd
                p = grid[i,j]
                p1 = grid[i+1,j]
                p2 = grid[i,j+1]
                p12 = calculate_mu_p12_pos(p, p1, p2, mu)
                grid[i+1,j+1] = p12
            else
                #i even, j odd
                p = grid[i,j+1]
                p1 = grid[i+1,j+1]
                p2 = grid[i,j]
                p12 = calculate_mu_p12_pos(p, p1, p2, mu)
                grid[i+1,j] = p12
            end
        else
            if i%2 == 1
                #i odd, j even
                p = grid[i+1,j]
                p1 = grid[i+1,j+1]
                p2 = grid[i,j]
                p12 = calculate_mu_p12_pos(p, p1, p2, mu)
                grid[i,j+1] = p12
            else
                #i even, j even
                p = grid[i,j]
                p1 = grid[i+1,j]
                p2 = grid[i,j+1]
                p12 = calculate_mu_p12_pos(p, p1, p2, mu)
                grid[i+1,j+1] = p12            
            end
        end
    end
end

on(p_move) do newval
    points_matrix[1,1] = p_move[]
    update_grid_cr_mu!(points_matrix,-1)
    points_matrix_isp = invStereographicProj.(points_matrix)
    points_matrix_isp_dual = discrete_dual(points_matrix_isp)

    obs_points_matrix[] = points_matrix
    obs_points_matrix_isp[] = points_matrix_isp
    obs_points_matrix_isp_dual[] = points_matrix_isp_dual
    
    update_edges!(obs_edges, points_matrix)
    update_edges!(obs_edges_isp, points_matrix_isp)
    update_edges!(obs_edges_isp_dual, points_matrix_isp_dual)
end


# draw quads
# for i in 1:n-1, j in 1:m-1
#     p, p1, p2, p12 = points_matrix[i,j],points_matrix[i+1,j],points_matrix[i,j+1],points_matrix[i+1,j+1]
#     quad_path = [p,p1,p12,p2,p]
#     lines!(ax,quad_path,color="grey")

    # s, s1, s2, s12 = points_matrix_isp[i,j],points_matrix_isp[i+1,j],points_matrix_isp[i,j+1],points_matrix_isp[i+1,j+1]
    # quad_face_3d = lift((P,P1,P12,P2) -> [P, P1, P12, P2], s, s1, s12 ,s2)
    # quad_path_3d = lift((P,P1,P12,P2) -> [P, P1, P12, P2, P], s, s1, s12, s2)
    # lines!(ax3d,quad_path_3d,color="grey")
    # poly!(ax3d, quad_face_3d, color=RGBA{Float32}(0, 0, 1, 0.1))

    # f, f1, f12, f2 = points_matrix_f_star[i,j],points_matrix_f_star[i+1,j],points_matrix_f_star[i+1,j+1],points_matrix_f_star[i,j+1]
    # quad_face_f = lift(f,f1,f12,f2) do P,P1,P12,P2
    #     [P,P1,P12,P2]
    # end
    # quad_path_f = lift(f,f1,f12,f2) do P,P1,P12,P2
    #     [P,P1,P12,P2,P]
    # end
    # lines!(ax3d,quad_path_f,color="grey")
    # poly!(ax3d,quad_face_f,color=RGBA{Float32}(0, 0, 1, 0.1)) 
    

    # circline is observable array of points that plot out a circle/line through p,p1,p2
#     circline = circle_or_line_coords_from_3pts(p,p1,p2)
#     lines!(ax,circline,color="green")
# end

is_dragging = Ref(false)
# A function to handle the logic when the mouse moves
on(events(fig).mouseposition) do pos
    # Only run the logic if the left mouse button is pressed and we are dragging
    if is_dragging[] && ispressed(fig, Mouse.left)
        # Convert the mouse position (in pixels) to the plot's data coordinates
        new_pos = Makie.mouseposition(ax.scene)
        
        # Update the Observable value, which automatically moves the point
        p_move[] = new_pos
        # println("Value in array: ", points_array[1,1][])
        # Consume the event to prevent the axis's default pan/zoom behavior
        return Consume(true) 
    end
    
    return Consume(false)
end

# A function to check if the drag *should* start
on(events(fig).mousebutton) do event
    # Only act on the left mouse button
    if event.button == Mouse.left
        # Start Drag: Check if the mouse was pressed AND the mouse is currently
        # close enough to the plot element (the red point).
        if event.action == Mouse.press
            # Use Makie's picking function to see what plot element is under the cursor
            picked_plot, picked_index = Makie.pick(ax.scene)
            
            # Check if the picked element is our scatter plot 'p_plot'
            if picked_plot == p_move_scat
                is_dragging[] = true
                return Consume(true)
            end
        
        # Stop Drag: If the mouse button is released, stop dragging.
        elseif event.action == Mouse.release
            is_dragging[] = false
            # We don't need to consume the release event
        end
    end
    return Consume(false)
end