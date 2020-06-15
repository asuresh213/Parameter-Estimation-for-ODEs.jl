using DifferentialEquations
using Plots


# restricted 3-body problem.

function prob_three_body!(du, u, p, t)
    x₁, y₁, x₂, y₂, m₁, m₂ = p
    r₁ = sqrt((u[1] - x₁)^2 + (u[2] - y₁)^2)
    r₂ = sqrt((u[1] - x₂)^2 + (u[2] - y₂)^2)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -(m₁/(r₁^3))*(u[1] - x₁) - (m₂/(r₂^3))*(u[1] - x₂)
    du[4] = -(m₁/(r₁^3))*(u[2] - y₁) - (m₂/(r₂^3))*(u[2] - y₂)
end

function ThreeBodyODE(u0, tspan :: NTuple{2, Float64}, p, dt_f=0.1)
    prob = ODEProblem(prob_three_body!, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    v = sol[1, :]
    w = sol[2, :]
    return sol.t, [v, w]
end

ThreeBodyODE(u0, tspan :: NTuple{2, Float64}, p; dt_f = 0.1) = ThreeBodyODE(u0, tspan, p, dt_f)
ThreeBodyODE(u0, tspan :: NTuple{2, Float64}; dt_f = 0.1) = ThreeBodyODE(u0, tspan, [40*rand()-20, 40*rand()-20, 40*rand()-20, 40*rand()-20, 100*rand(), 100*rand()], dt_f)
ThreeBodyODE(tspan :: NTuple{2, Float64}, p; dt_f = 0.1) = ThreeBodyODE(20 .* rand(4) .-10, tspan, p, dt_f)





# example function call
u0 = [1.0, 0.0, 0.0, -5.00158510637908252240537862224]
tspan = (0.0, 55.0652165601579625588917206249)
p = [5.0, 0.0, -5.0, 0.0, 100.0, 100.0]
testtime, testsol = ThreeBodyODE(u0, tspan, p; dt_f=0.01)
scatter(testsol[1], testsol[2])


# plotting and making a nice gif of the plot!
plot(testsol[1], testsol[2])
n = 20
@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 2, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

@gif for i in 1:length(testtime)
    circleplot(testsol[1], testsol[2], i, line_z = 1:n, cbar = false, c = :reds, framestyle = :none)
end when i>(20) && mod1(i,10) == 5
