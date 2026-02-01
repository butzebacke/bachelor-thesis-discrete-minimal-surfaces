using GLMakie
using Observables
using LinearAlgebra
using Colors
using GeometryBasics


function update_g_cr_evol(xs,ys)
    n,m = length(xs), length(ys)
    for i in 1:n
        g[(i,1)] = xs[i]
    end
    for j in 2:m
        g[(1,j)] = ys[j]
    end
    for j in 2:m
        for i in 2:n
            p = g[(i-1,j-1)]
            p1 = g[(i,j-1)]
            p2 = g[(i-1,j)]
            p12 = calculate_p12_conformal_square(p,p1,p2)
            g[(i,j)] = p12
        end
    end
end

function get_weierstrass_edge(g,idx1,idx2)
    gstart = g[idx1]
    gend = g[idx2]

    edge = 1/2 * real.([1 - gstart*gend, im*(1 + gstart*gend), gstart + gend] ./ (gend - gstart))

    if (idx1[2] == idx2[2]) # \Delta_1 edge
        return edge
    else # \Delta_2 edge
        return -1 .* edge
    end
end

function weierstrass_construct_f(g)
    f[(1,1)] = Point3f(0)
    for j in 2:m
        fstart = f[(1,j-1)]
        del2f = get_weierstrass_edge(g,(1,j-1),(1,j))
        f[(1,j)] = fstart + del2f
    end
    for i in 2:n, j in 1:m
        fstart = f[(i-1,j)]
        del1f = get_weierstrass_edge(g,(i-1,j),(i,j))
        f[(i,j)] = fstart + del1f
    end

    # balance
    shift = Point3f(0)
    k = 0
    for v in values(f) 
        shift += v
        k += 1
    end
    shift /= k
    for idx in keys(f) f[idx] -= shift end
end

function complex_to_p2f(z::ComplexF64)
    return Point2f((real(z),imag(z)))
end

function p2f_to_complex(p::Point2f)
    return p[1] + p[2] * im
end

function calculate_p12_conformal_square(p::ComplexF64, p1::ComplexF64, p2::ComplexF64) 
    return (p2*(p-p1) - p1*(p2-p)) / ((2*p - p1 - p2))
end

function plot_g(ax; vertexcolor=:black, edgecolor=:blue)
    # plot edges
    for i in 1:n
        segs = [complex_to_p2f(g[(i,j)]) for j in 1:m]
        lines!(ax, segs, color=edgecolor)
    end
    for j in 1:m 
        segs = [complex_to_p2f(g[(i,j)]) for i in 1:n]
        lines!(ax, segs, color=edgecolor)
    end
    scatter!(ax, complex_to_p2f.(values(g)), color=vertexcolor, markersize=6)
end

function plot_f(ax; showfaces=true, vertexcolor=:black, edgecolor=:black, facecolor=:steelblue, vertexsize=8, edgesize=2.5)
    if showfaces
        for i in 1:n-1,j in 1:m-1
            face = [f[(i,j)],f[(i+1,j)],f[(i+1,j+1)],f[(i,j+1)]]
            poly!(ax, face, color=(facecolor, 0.5), transparency=true)
        end
    end
    for i in 1:n
        segs = [f[(i,j)] for j in 1:m]
        lines!(ax, segs, color=(edgecolor, 0.8), linewidth=edgesize)
    end
    for j in 1:m 
        segs = [f[(i,j)] for i in 1:n]
        lines!(ax, segs, color=(edgecolor, 0.8), linewidth=edgesize)
    end
    scatter!(ax, collect(values(f)), color=vertexcolor, markersize=vertexsize)
end

# formulas taken from Bobenko.1996
function discrete_exp1996_init_vals(n,m,N)
    alpha = 2*pi/N
    x = Vector{ComplexF64}()
    y = Vector{ComplexF64}()

    for i in 1:n
        rho = 2*asinh(sin(alpha/2))
        xi = exp(rho * (i-1))
        push!(x,xi)
    end
    for j in 1:m
        yj = exp(im*alpha*(j-1))
        push!(y,yj)
    end

    return x,y
end

function id_init_vals(n,m,minx,miny,distx,disty)
    x = Vector{ComplexF64}()
    y = Vector{ComplexF64}()
    minx,miny = Float64(minx),Float64(miny)
    distx,disty = Float64(distx),Float64(disty)
    z0 = minx + miny*im
    for i in 1:n
        push!(x,z0 + (i-1)*distx)
    end
    for j in 1:m
        push!(y,z0 + im*(j-1)*disty)
    end
    return x,y
end

function nice_catenoid_init_vals()
    x,y = discrete_exp1996_init_vals(13, 17, 16)
    return x*0.1,y*0.1
end

x_init_vals, y_init_vals = id_init_vals(3,3,-0.5,-0.5,0.5,0.5)
x_init_vals[2] -= 0.2
y_init_vals[2] += 0.2
n = length(x_init_vals)
m = length(y_init_vals)

g = Dict{Tuple{Int,Int},ComplexF64}()
update_g_cr_evol(x_init_vals,y_init_vals)

f = Dict{Tuple{Int,Int},Point3f}()
weierstrass_construct_f(g)

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

plot_g(ax; edgecolor=:grey)
plot_f(ax3d)

display(fig)
