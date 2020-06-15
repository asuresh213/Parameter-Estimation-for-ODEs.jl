using DifferentialEquations
using Plots

prob_lotka_volterra! = @ode_def begin
    dx = a*x - b*x*y
    dy = -c*y + d*x*y
end a b c d

function LotkaVolterraODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}, dt_f :: Number)
    prob = ODEProblem(prob_lotka_volterra!, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    prey = sol[1, :]
    pred = sol[2, :]
    return sol.t, [prey, pred]
end

# given u0, tspan, p, dt_f
LotkaVolterraODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}, dt_f :: Number) = LotkaVolterraODE(u0, tspan, [i for i in p], dt_f)
LotkaVolterraODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 0.1) = LotkaVolterraODE(u0, tspan, p, dt_f)
LotkaVolterraODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}; dt_f = 0.1) = LotkaVolterraODE(u0, tspan, p, dt_f)

# given tspan
# randomized IC with pred and prey between 0, 3.
# randomzed parameters, each component between 0 and 5.
LotkaVolterraODE(tspan :: NTuple{2, Number}; dt_f = 0.1) = LotkaVolterraODE(3*rand(2) .+ 0.5, tspan, 5*rand(4) .+ 0.2, dt_f)

# given u0, tspan
# randomzed parameters, each component between 0 and 5.
LotkaVolterraODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}; dt_f = 0.1) = LotkaVolterraODE(u0, tspan, 5*rand(4) .+ 0.2, dt_f)

# given tspan, p
# randomzed u0, each component between 0 and 5.
LotkaVolterraODE(tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 0.1) = LotkaVolterraODE(3*rand(2) .+ 0.5, tspan, p, dt_f)

function LotkaVolterraSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}, dt_f = 0.1, noise_factor = 0.15)
    function lotka_volterra_noise!(du, u, p, t)
        du[1] = noise_factor*u[1]
        du[2] = noise_factor*u[2]
    end
    prob = SDEProblem(prob_lotka_volterra!, lotka_volterra_noise!, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    prey = sol[1, :]
    pred = sol[2, :]
    return sol.t, [prey, pred]
end

# given u0, tspan, p, dt_f
LotkaVolterraSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}, dt_f :: Number, noise_factor :: Number) = LotkaVolterraSDE(u0, tspan, [i for i in p], dt_f, noise_factor)
LotkaVolterraSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 0.1, noise_factor = 0.15) = LotkaVolterraSDE(u0, tspan, p, dt_f, noise_factor)
LotkaVolterraSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}; dt_f = 0.1, noise_factor = 0.15) = LotkaVolterraSDE(u0, tspan, p, dt_f, noise_factor)

# given tspan
# randomized IC with pred and prey between 0, 3.
# randomzed parameters, each component between 0 and 5.
LotkaVolterraSDE(tspan :: NTuple{2, Number}; dt_f = 0.1, noise_factor = 0.15) = LotkaVolterraSDE(3*rand(2) .+ 0.5, tspan, 5*rand(4) .+ 0.2, dt_f, noise_factor)

# given u0, tspan
# randomzed parameters, each component between 0 and 5.
LotkaVolterraSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}; dt_f = 0.1, noise_factor = 0.15) = LotkaVolterraSDE(u0, tspan, 5*rand(4) .+ 0.2, dt_f, noise_factor)

# given tspan, p
# randomzed u0, each component between 0 and 5.
LotkaVolterraSDE(tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 0.1, noise_factor = 0.15) = LotkaVolterraSDE(3*rand(2) .+ 0.5, tspan, p, dt_f, noise_factor)


# example function call

testtime, testsol = LotkaVolterraODE([1.0, 1.0], (0.0, 10.0), (2.0, 1.0, 1.0, 2.0); dt_f = 0.01)
scatter(testtime, testsol[1])
scatter!(testtime, testsol[2])
