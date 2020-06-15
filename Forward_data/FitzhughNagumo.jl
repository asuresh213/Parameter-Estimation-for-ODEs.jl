using DifferentialEquations
using Plots


prob_fitz_nagumo! = @ode_def begin
    dv = v - (1/3)*v^3 - w + (1/10)*(5 + sin((3.14*t)/10))
    dw = (v + a - b*w)/τ
end a b τ

function FitzhughNagumoODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}, dt_f :: Number)
    prob = ODEProblem(prob_fitz_nagumo!, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    v = sol[1, :]
    w = sol[2, :]
    return sol.t, [v, w]
end

# given u0, tspan, p, dt_f
FitzhughNagumoODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}, dt_f :: Number) = FitzhughNagumoODE(u0, tspan, [i for i in p], dt_f)
FitzhughNagumoODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 1.0) = FitzhughNagumoODE(u0, tspan, p, dt_f)
FitzhughNagumoODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}; dt_f = 1.0) = FitzhughNagumoODE(u0, tspan, p, dt_f)

# given tspan
# randomized IC with pred and prey between 0, 3.
# randomzed parameters, each component between 0 and 5.
FitzhughNagumoODE(tspan :: NTuple{2, Number}; dt_f = 1.0) = FitzhughNagumoODE(3*rand(2) .+ 0.5, tspan, [2*rand(), 3*rand(), 30*rand()] .+ 0.2, dt_f)

# given u0, tspan
# randomzed parameters, each component between 0 and 5.
FitzhughNagumoODE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}; dt_f = 1.0) = FitzhughNagumoODE(u0, tspan, [2*rand(), 3*rand(), 30*rand()] .+ 0.2, dt_f)

# given tspan, p
# randomzed u0, each component between 0 and 5.
FitzhughNagumoODE(tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 1.0) = FitzhughNagumoODE(3*rand(2) .+ 0.5, tspan, p, dt_f)

function FitzhughNagumoSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}, dt_f = 1.0, noise_factor = 0.05)
    function fitz_nagumo_noise!(du, u, p, t)
        du[1] = noise_factor*u[1]
        du[2] = noise_factor*u[2]
    end
    prob = SDEProblem(prob_fitz_nagumo!, fitz_nagumo_noise!, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    prey = sol[1, :]
    pred = sol[2, :]
    return sol.t, [prey, pred]
end

# given u0, tspan, p, dt_f
FitzhughNagumoSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}, dt_f :: Number, noise_factor :: Number) = FitzhughNagumoSDE(u0, tspan, [i for i in p], dt_f, noise_factor)
FitzhughNagumoSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 1.0, noise_factor = 0.05) = FitzhughNagumoSDE(u0, tspan, p, dt_f, noise_factor)
FitzhughNagumoSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}, p :: Vector{Float64}; dt_f = 1.0, noise_factor = 0.05) = FitzhughNagumoSDE(u0, tspan, p, dt_f, noise_factor)

# given tspan
# randomized IC with pred and prey between 0, 3.
# randomzed parameters, each component between 0 and 5.
FitzhughNagumoSDE(tspan :: NTuple{2, Number}; dt_f = 1.0, noise_factor = 0.05) = FitzhughNagumoSDE(3*rand(2) .+ 0.5, tspan, [2*rand(), 3*rand(), 30*rand()] .+ 0.2, dt_f, noise_factor)

# given u0, tspan
# randomzed parameters, each component between 0 and 5.
FitzhughNagumoSDE(u0 :: Vector{Float64}, tspan :: NTuple{2, Number}; dt_f = 1.0, noise_factor = 0.05) = FitzhughNagumoSDE(u0, tspan, [2*rand(), 3*rand(), 30*rand()] .+ 0.2, dt_f, noise_factor)

# given tspan, p
# randomzed u0, each component between 0 and 5.
FitzhughNagumoSDE(tspan :: NTuple{2, Number}, p :: NTuple{4, Number}; dt_f = 1.0, noise_factor = 0.05) = FitzhughNagumoSDE(3*rand(2) .+ 0.5, tspan, p, dt_f, noise_factor)



# example function call
testtime, testsol = FitzhughNagumoODE([1.0, 1.0], (0.0, 100.0), (0.7, 1.1, 14.0))
scatter(testtime, testsol[1])
scatter!(testtime, testsol[2])
