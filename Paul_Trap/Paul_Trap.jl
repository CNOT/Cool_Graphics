using OrdinaryDiffEq
using GLMakie

#Parameters, we can later replace them with sliders
Ω = 3.0
aᵤ = 0.1
qᵤ = 0.8
p = [Ω,aᵤ,qᵤ]

#Initial conditions
u₀ = 2*rand(6)-ones(6)
tspan = (0.0, 10π)

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

xs = [sol.u[i][1] for i in 1:length(sol.u)]
ys = [sol.u[i][2] for i in 1:length(sol.u)]
zs = [sol.u[i][3] for i in 1:length(sol.u)]
vxs = [sol.u[i][4] for i in 1:length(sol.u)]
vys = [sol.u[i][5] for i in 1:length(sol.u)]
vzs = [sol.u[i][6] for i in 1:length(sol.u)]


# Plot
# using Plots
# plot(sol, vars=(1,2) , linewidth=2, title="Paul Trap", xaxis="x", yaxis="y")

# Animation

points = Observable(Point2f[(u₀[1],u₀[2])])
colors = Observable(Float64[])
frames = 1:length(xs)
max_x = max(abs.(xs)...) 
max_y = max(abs.(ys)...)
max_z = max(abs.(zs)...)
fig, ax = scatter(points, color = :blue, markersize = [1])
limits!(ax, -1.2*max_x, 1.2*max_x, -1.2*max_y, 1.2*max_y)
record(fig, "Paul_Trap/PaulTrap_animation.mp4", frames;
framerate=30) do frame
    new_point = Point2f(xs[frame], ys[frame])
    points[] = push!(points[], new_point)
    colors = [(:blue, (2^i)/2^(frame+1)) for i in 1:frame+1]
    scatterlines!(points, color=:white, markersize=repeat([2], frame + 1))
    scatterlines!(points, color = colors, markersize = repeat([1],frame+1))

    # arrows!([xs[frame]],[ys[frame]],[vxs[frame]],[vys[frame]],lengthscale = 0.2, arrowcolor = colors)
end