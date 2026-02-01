using GLMakie
using Observables
using LinearAlgebra
using Colors
using GeometryBasics

function complex_to_p2f(z::ComplexF64)
    return Point2f((real(z),imag(z)))
end

function p2f_to_complex(p::Point2f)
    return p[1] + p[2] * im
end

function orientation_is_positive(z1,z2,z3)
    return imag((z3 - z1) * conj(z2 - z1)) > 0
end

function calculate_r2_orthogonal_neighbor(m1,m2,r1)
    d = abs(m2-m1)
    r2sq = d^2 - r1^2
    r2sq > 0 ? r2 = sqrt(r2sq) : error("circle with midpoints: ",m1 ," , ", m2, " are inside each other")
    return r2
end

function calculate_intersections(c1::Tuple{Int,Int}, c2::Tuple{Int,Int})
    m1, m2 = g_w[c1], g_w[c2]
    r1, r2 = rs[c1], rs[c2]        
    d = abs(m1 - m2)
    u = (m2 - m1) / d # unit vector from m1 to m2
    a = (r1^2 - r2^2 + d^2) / (2*d) 
    h = sqrt(Complex(r1^2 - a^2))
    z1 = m1 + a*u + im*h*u
    z2 = m1 + a*u - im*h*u
    if orientation_is_positive(m1,m2,z1)
        return (z1,z2)
    else 
        return (z2,z1)
    end
end

function set_circle(c::Tuple{Int,Int},m,r)
    i,j = c[1],c[2]
    if (i + j) % 2 == 1
        throw(DomainError(c))
    end
    if r >= 0
        rs[c] = r 
        g_w[c] = m
    else
        throw(DomainError(r))
    end
    # Adding intersecions
    if haskey(g_w, c.+(1,1))
        (z1,z2) = calculate_intersections(c,c.+(1,1))
        g_b[(i,j+1)] = z1
        g_b[(i+1,j)] = z2
    end
    if haskey(g_w, c.+(-1,1))
        (z1,z2) = calculate_intersections(c,c.+(-1,1))
        g_b[(i-1,j)] = z1
        g_b[(i,j+1)] = z2
    end
    if haskey(g_w, c.+(-1,-1))
        (z1,z2) = calculate_intersections(c,c.+(-1,-1))
        g_b[(i,j-1)] = z1
        g_b[(i-1,j)] = z2
    end
    if haskey(g_w, c.+(1,-1))
        (z1,z2) = calculate_intersections(c,c.+(1,-1))
        g_b[(i+1,j)] = z1
        g_b[(i,j-1)] = z2
    end
end

# calculate midpoint of c that is tangent to c1 and c2 in the intersections with cmid and orientation m1,m2,m positive
function set_tangent_c(c1,c2)
    if abs(c1[1] - c2[1]) == 2 && abs(c1[2] - c2[2]) == 2   
        m1, m2 = g_w[c1], g_w[c2]
        v = Int.((c2 .- c1) ./ 2)
        cmid = c1 .+ v
        # intersections s.t. b1 = c1 \cap cmid and b2 = c2 \cap cmid, avoid matrix calculation for speed I think
        if (v[1] * v[2] > 0)
            b1 = cmid .+ (-v[1],0)
            b2 = cmid .+ (0,v[2])
        else
            b1 = cmid .+ (0,-v[2])
            b2 = cmid .+ (v[1],0)
        end
        cnew = b1 .+ (b1 .- c1)
        z1, z2 = g_b[b1], g_b[b2]
        d1 = z1-m1
        d2 = z2-m2
        d1real = [real(d1),imag(d1)]
        d2real = [real(d2),imag(d2)]
        # solve lgs for intersection of m1 + t*d1 = m2 + s*d2
        A = hcat(d1real,-d2real)
        b = [real(m2-m1),imag(m2-m1)]
        ts = A \ b
        mnew = m1 + ts[1]*d1
        rnew = abs(mnew - z1)
        set_circle(cnew,mnew,rnew)
    else
        error("c1 and c2 dont have a common tangent circle")
    end
end

function update_g_diag_evol(ms,r0)
    n = length(ms)
    set_circle((1,1),ms[1],r0)
    # set diag circles
    for i in 2:n 
        m1 = g_w[(i-1,i-1)]
        m2 = ms[i]
        r1 = rs[(i-1,i-1)]
        r2 = calculate_r2_orthogonal_neighbor(m1,m2,r1)
        set_circle((i,i),m2,r2)
    end

    for k in 0:2:n-1
        c_up_start = (1, 1 + k)
        c_low_start = (1 + k, 1)
        while c_up_start[2] <= n-2
            set_tangent_c(c_up_start, c_up_start .+ (2,2))
            set_tangent_c(c_low_start .+ (2,2), c_low_start)
            c_up_start = c_up_start .+ (1,1)
            c_low_start = c_low_start .+ (1,1)
        end
    end
end

function circle_center_3d(p1, p2, p3; p4=nothing, atol=1e-5)
    v2 = p2 - p1
    v3 = p3 - p1

    n = cross(v2, v3)
    nn = norm(n)
    nn == 0 && error("p1,p2,p3 are collinear (no unique circle).")

    e1 = v2 / norm(v2)      # in-plane unit vector
    e3 = n / nn             # plane normal unit vector
    e2 = cross(e3, e1)      # in-plane unit vector orthogonal to e1

    x2, y2 = dot(v2, e1), dot(v2, e2)
    x3, y3 = dot(v3, e1), dot(v3, e2)
    d = 2 * (x2*y3 - y2*x3)

    abs(d) < atol && error("Nearly collinear / numerically unstable.")

    r2 = x2^2 + y2^2
    r3 = x3^2 + y3^2

    cx = (r2*y3 - r3*y2) / d
    cy = (r3*x2 - r2*x3) / d

    center = p1 + cx*e1 + cy*e2

    # Optional consistency check with p4
    if p4 !== nothing
        R = norm(center - p1)
        err_radius = abs(norm(center - p4) - R)
        err_plane  = abs(dot(p4 - p1, e3))
        (err_radius > atol || err_plane > atol) && @warn(
            "p4 not consistent with same circle/plane within tolerance",
            err_radius=err_radius, err_plane=err_plane
        )
    end

    return center
end
# get radius of sphere at c
function get_weierstrass_radius(c::Tuple{Int,Int})
    # black vertex for calculation
    b = haskey(g_b,c.+(1,0)) ? c.+(1,0) : c.-(1,0)
    center_val = g_w[c]
    intersection_val = g_b[b]
    r = abs((1 + abs2(center_val) - abs2(center_val - intersection_val)) /
    (2*abs(center_val - intersection_val)))
    return r
end 

function get_weierstrass_s_center_difference_vector(c1::Tuple{Int,Int},c2::Tuple{Int,Int})
    sgn = (c2[1] == c1[1]) ? -1 : 1
    b = Int.((c1 .+ c2) ./ 2)
    b_val = g_b[b]
    c1_val = g_w[c1]
    c2_val = g_w[c2]
    A = (rf[c1] + rf[c2]) / (1 + abs2(b_val))
    B = (conj(c2_val) - conj(c1_val)) / abs(c2_val - c1_val)
    Phi = [1 - b_val^2, im*(1 + b_val^2), 2*b_val]
    return Point3f((sgn .* real.(A * B .* Phi))...)
end

# first calculate all s-vertices via weierstrass and their radii, then fill intersections, then calculate c-vertices and their radii
# for now this only works if the left n_min and n_max are both s-vertices
function calculate_weierstrass_smin_from_square()
    
    # fill radii
    for i in 1:2:n, j in 1:2:n
        rf[(i,j)] = get_weierstrass_radius((i,j))
    end

    f_w[(1,1)] = Point3f(0)
    # fill left column of s-centers
    for j in 3:2:n
        c1 = (1,j-2)
        c2 = (1,j)
        del_f = get_weierstrass_s_center_difference_vector(c1,c2)
        f_w[c2] = f_w[c1] .+ del_f
    end

    for i in 3:2:n, j in 1:2:n
        c1 = (i-2,j)
        c2 = (i,j)
        del_f = get_weierstrass_s_center_difference_vector(c1,c2)
        f_w[c2] = f_w[c1] .+ del_f
    end

    for i in 1:2:n, j in 2:2:n-1
        c1 = (i,j-1)
        c2 = (i,j+1)
        b = (i,j)
        v = f_w[c2] .- f_w[c1]
        v_norm = v ./ norm(v)
        r1 = rf[c1]
        f_b[b] = f_w[c1] .+ r1.*v_norm
    end
    for j in 1:2:n, i in 2:2:n-1
        c1 = (i-1,j)
        c2 = (i+1,j)
        b = (i,j)
        v = f_w[c2] .- f_w[c1]
        v_norm = v ./ norm(v)
        r1 = rf[c1]
        f_b[b] = f_w[c1] .+ r1.*v_norm
    end

    # c-centers 
    for i in 2:2:n-1, j in 2:2:n-1
        p1 = f_b[(i+1,j)]
        p2 = f_b[(i,j+1)]
        p3 = f_b[(i-1,j)]
        center = circle_center_3d(p1,p2,p3)
        f_w[(i,j)] = center
        rf[(i,j)] = norm(center .- p1)
    end
end

function plot_g(ax; show_intersections = false, show_kites = false, n_segs = 50, kitecolor=:darkgoldenrod, circlecolor=:steelblue, edgewidth=2.2, vertexsize=6, nmin=1, nmax=n, mmin=1, mmax=n)
    phi = range(0,2*pi; length=n_segs)
    midpoints = Vector{Point2f}()
    for k in keys(g_w)
        if (nmin <= k[1] <= nmax && mmin <= k[2] <= mmax) 
            m = complex_to_p2f(g_w[k])
            r = rs[k]
            lines!(ax, m[1] .+ r*cos.(phi), m[2] .+ r*sin.(phi); color=circlecolor, linewidth = edgewidth)
            push!(midpoints,m)
            if show_kites
                if haskey(g_b, k.+(1,0))
                    b = complex_to_p2f(g_b[(k.+(1,0))])
                    lines!(ax, [m,b]; color=kitecolor, linewidth = edgewidth)
                end
                if haskey(g_b, k.+(0,1))
                    b = complex_to_p2f(g_b[(k.+(0,1))])
                    lines!(ax, [m,b]; color=kitecolor, linewidth = edgewidth)
                end
                if haskey(g_b, k.+(-1,0))
                    b = complex_to_p2f(g_b[(k.+(-1,0))])
                    lines!(ax, [m,b]; color=kitecolor, linewidth = edgewidth)
                end
                if haskey(g_b, k.+(0,-1))
                    b = complex_to_p2f(g_b[(k.+(0,-1))])
                    lines!(ax, [m,b]; color=kitecolor, linewidth = edgewidth)
                end
            end
        end
    end
    # scatter!(ax, midpoints; markersize=vertexsize, color=:white, strokecolor=:black, strokewidth=1)
    if show_intersections
        intersections = collect(complex_to_p2f.(values(g_b)))
        scatter!(ax, intersections; markersize=vertexsize, color=:black)
    end
end

function plot_f(ax; vertexcolor=:black, spherecolor=:steelblue, spherealpha=0.8, circlecolor=:darkgoldenrod, circlealpha=0.8, N=30, showspheres=true, showcircles=true, nmin=1, mmin=1, nmax=n, mmax=n)
    for (c, m) in f_w 
        if (nmin <= c[1] <= nmax) && (mmin <= c[2] <= mmax)
            if isodd(c[1]) && isodd(c[2]) && showspheres
                haskey(rf, c) || continue
                r = rf[c]
                # skip degenerate / invalid spheres
                (isfinite(r) && r > 1e-6) || continue

                # make sure we have Float32 points/radii (Makie likes that)
                center = Point3f(Float32.(m)...)
                rad    = Float32(r)

                mesh!(ax3d, Sphere(center, rad), color= (spherecolor, spherealpha), transparency=true)
                # optional: wireframe outline
                # wireframe!(ax, Sphere(center, rad), color=:black, linewidth=0.5)
            elseif iseven(c[1]) && iseven(c[2]) && showcircles
                haskey(rf, c) || continue
                r = rf[c]

                (isfinite(r) && r > 1e-6) || continue

                if haskey(f_b, c.+(1,0)) 
                    b1 = f_b[c.+(1,0)]
                else 
                    b1 = f_b[c.+(-1,0)]
                end
                if haskey(f_b, c.+(0,1))
                    b2 = f_b[c.+(0,1)]
                else
                    b2 = f_b[c.+(0,-1)]
                end

                # build ONB u,v,n 
                u = b1 - m
                u /= norm(u)

                ng = cross(u, b2-m)
                ng /= norm(ng)
                v = cross(ng,u)

                theta = range(0,2*pi; length = N)
                pts = [m + r*(cos(t)*u + sin(t)*v) for t in theta]
                poly!(ax, pts; color = (circlecolor, circlealpha))
            end
        end
    end
end

function id_init_m(n,xmin,ymin,r0)
    m = Vector{ComplexF64}()
    z = xmin + im*ymin
    for i in 1:n
        push!(m,z)
        z += r0*(1+1*im)
    end
    return m
end

function exp_init_m(N)
    rho = 2*pi / (N-1) 
    alpha = atanh(1/2 * abs(1 - exp(2*im*rho)))
    m = Vector{ComplexF64}()
    push!(m, 1.0+0.0*im)
    for i in 1:N-1
        z = exp((alpha + im*rho)*i)
        push!(m,z)
    end
    r0 = sin(rho)
    return m, r0
end

fig = Figure()
ax = Axis(fig[1,1], aspect=DataAspect())

isaxlabelvisible = false
ax3d = Axis3(fig[1,2], xlabelvisible = isaxlabelvisible,
    ylabelvisible = isaxlabelvisible,
    zlabelvisible = isaxlabelvisible,
    xticksvisible = isaxlabelvisible,
    yticksvisible = isaxlabelvisible,
    zticksvisible = isaxlabelvisible,
    xticklabelsvisible = isaxlabelvisible,
    yticklabelsvisible = isaxlabelvisible,
    zticklabelsvisible = isaxlabelvisible)
ax3d.aspect = :data

g_w = Dict{Tuple{Int,Int},ComplexF64}()
g_b = Dict{Tuple{Int,Int},ComplexF64}()
rs = Dict{Tuple{Int,Int},Float64}()

f_w = Dict{Tuple{Int,Int},Point3f}()
f_b = Dict{Tuple{Int,Int},Point3f}()
rf = Dict{Tuple{Int,Int},Float64}()

r0 = 0.05
init_m = id_init_m(21,-1,-1,r0)
n = length(init_m)

update_g_diag_evol(init_m, r0)
calculate_weierstrass_smin_from_square()
plot_g(ax; show_intersections=false, show_kites=false, vertexsize=8)
plot_f(ax3d; spherecolor=:firebrick1, spherealpha=0.5, circlecolor=:steelblue, circlealpha=1.0, showspheres=true,showcircles=true)
println(f_w)

display(fig)