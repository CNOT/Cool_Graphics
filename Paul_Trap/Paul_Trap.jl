using OrdinaryDiffEq
using GLMakie
using CairoMakie
#Parameters, we can later replace them with sliders
Ω = 3.0
aᵤ = 0.1
qᵤ = 0.8
p = [Ω,aᵤ,qᵤ]

#Initial conditions
u₀ = 2*rand(6)-ones(6)
tspan = (0.0, 20π)

# Trap function
function PaulTrap!(du,u,p,t)
    for i in 1:3
        du[i] = u[i+3]
        du[i+3] = (Ω^2)/4 *  (2qᵤ*cos(Ω*t)-aᵤ)*u[i]
    end
end

# ODE problem and solution
prob = ODEProblem(PaulTrap!, u₀, tspan)
sol = solve(prob, Tsit5())

xs = [sol.(0:0.1:20π)[i][1] for i in 1:length(sol.(0:0.1:20π))]
ys = [sol.(0:0.1:20π)[i][2] for i in 1:length(sol.(0:0.1:20π))]
zs = [sol.(0:0.1:20π)[i][3] for i in 1:length(sol.(0:0.1:20π))]
vxs = [sol.(0:0.1:20π)[i][4] for i in 1:length(sol.(0:0.1:20π))]
vys = [sol.(0:0.1:20π)[i][5] for i in 1:length(sol.(0:0.1:20π))]
vzs = [sol.(0:0.1:20π)[i][6] for i in 1:length(sol.(0:0.1:20π))]


# Plot

plt_xy , ax_fig_xy = scatterlines(xs,ys, markersize = 1,color = [vxs[i]^2 + vys[i]^2 for i in 1:length(xs)], colormap = :imola)
plt_yz, ax_fig_yz = scatterlines(ys, zs, markersize=1, color=[vzs[i]^2 + vys[i]^2 for i in 1:length(xs)], colormap=:imola)
plt_zx, ax_fig_zx = scatterlines(zs, xs, markersize=1, color=[vxs[i]^2 + vzs[i]^2 for i in 1:length(xs)], colormap=:imola)


# Animation

frames = 16:length(xs)
max_x = max(abs.(xs)...) 
max_y = max(abs.(ys)...)
max_z = max(abs.(zs)...)

# XY plane
points = Observable(Point2f[(xs[i],ys[i]) for i in 1:15])
fig_xy, ax_xy = scatterlines(points, color=15:-1:1, colormap=:viridis, markersize=0.5*(1:15))
limits!(ax_xy, -1.2*max_x, 1.2*max_x, -1.2*max_y, 1.2*max_y)
record(fig_xy, "Paul_Trap/PaulTrap_xy_animation.mp4", frames;
framerate=30) do frame
    new_point = Point2f(xs[frame], ys[frame])
    popfirst!(points[])
    points[] = push!(points[], new_point)
    scatterlines!(points, color=15:-1:1, colormap=:viridis, markersize=0.5*(1:15))

    # arrows!([xs[frame]],[ys[frame]],[vxs[frame]],[vys[frame]],lengthscale = 0.2, arrowcolor = colors)
end


# YZ Plane
points = Observable(Point2f[(ys[i], zs[i]) for i in 1:15])
fig_yz, ax_yz = scatterlines(points, color=15:-1:1, colormap=:viridis, markersize=0.5*(1:15))
limits!(ax_yz, -1.2 * max_y, 1.2 * max_y, -1.2 * max_z, 1.2 * max_z)
record(fig_yz, "Paul_Trap/PaulTrap_yz_animation.mp4", frames;
    framerate=30) do frame
    new_point = Point2f(ys[frame], zs[frame])
    popfirst!(points[])
    points[] = push!(points[], new_point)
    scatterlines!(points, color=15:-1:1, colormap=:viridis, markersize=0.5*(1:15))

    # arrows!([xs[frame]],[ys[frame]],[vxs[frame]],[vys[frame]],lengthscale = 0.2, arrowcolor = colors)
end

# ZX Plane
points = Observable(Point2f[(zs[i], xs[i]) for i in 1:15])
fig_zx, ax_zx = scatterlines(points, color=15:-1:1, colormap=:viridis, markersize=0.5*(1:15))
limits!(ax_zx, -1.2 * max_z, 1.2 * max_z, -1.2 * max_x, 1.2 * max_x)
record(fig_zx, "Paul_Trap/PaulTrap_zx_animation.mp4", frames;
    framerate=30) do frame
    new_point = Point2f(zs[frame], xs[frame])
    popfirst!(points[])
    points[] = push!(points[], new_point)
    scatterlines!(points, color=15:-1:1, colormap=:viridis, markersize=0.5*(1:15))

    # arrows!([xs[frame]],[ys[frame]],[vxs[frame]],[vys[frame]],lengthscale = 0.2, arrowcolor = colors)
end

