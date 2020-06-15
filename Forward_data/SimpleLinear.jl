using DifferentialEquations
using Plots

# Solving a simlple ode system
        # y' = α*y
        # with
            # u0 = Initial condition (number)
            # tspan = Time Span of integration (tuple)
            # p = Parameter α (number)
            # dt_f = discrete time step value


# ODE Objective function
prob_ode_linear!(u, p, t) = p*u


# Solving the ODE problem using DifferentialEquations.jl with user given u0, tspan, p, dt_f
function SimpleLinearODEModel(u0 :: Number, tspan :: Tuple{Number, Number}, p :: Number, dt_f :: Number)
    prob = ODEProblem(prob_ode_linear!, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    return sol.t, sol.u
end

# Solving the ODE problem using DifferentialEquations.jl with user given tspan
# u0, p are assigned (uniform) random values ∈ [-2.5,2.5]
# dt_f is set to 0.1 as a default measure.

SimpleLinearODEModel(tspan :: Tuple{Number, Number}; dt_f = 0.1) = SimpleLinearODEModel(5*rand() - 2.5, tspan, 5*rand() - 2.5, dt_f)

# Solving the ODE problem using DifferentialEquations.jl with user given u0, tspan
# p is assigned a (uniform) random value ∈ [-2.5,2.5]
# dt_f is set to 0.1 as a default measure.

SimpleLinearODEModel(u0 :: Number, tspan :: Tuple{Number, Number}; dt_f = 0.1) = SimpleLinearODEModel(u0, tspan, typeof(u0)(5*rand()-2.5), typeof(u0)(dt_f))

# Solving the ODE problem using DifferentialEquations.jl with user given tspan, p
# u0 is assigned a (uniform) random value ∈ [-2.5,2.5]
# dt_f is set to 0.1 as a default measure.

SimpleLinearODEModel(tspan :: Tuple{Number, Number}, p :: Number; dt_f = 0.1) = SimpleLinearODEModel(typeof(p)(5*rand() - 2.5), tspan, p, typeof(p)(dt_f))


# Solving the ODE problem using DifferentialEquations.jl with user given u0, tspan, p
# dt_f is set to 0.1 as a default measure.

SimpleLinearODEModel(u0 :: Number, tspan :: Tuple{Number, Number}, p :: Number; dt_f = 0.1) = SimpleLinearODEModel(u0, tspan, p, typeof(p)(dt_f))

# ------------------------------------------------------------------------------

# Solving the linear SDE of the form
    # y' = f + gdW = α*y + noise_factor*(y)dW
    # with
        # u0 = Initial condition (number)
        # tspan = Time Span of integration (tuple)
        # p = Parameter α (number)
        # dt_f = discrete time step value for f
        # noise_factor = factor by which the noise gets amplified

# Solving the SDE problem using DifferentialEquations.jl with user given u0, tspan, p, dt_f, noise_factor
function SimpleLinearSDEModel(u0 :: Number, tspan :: Tuple{Number, Number}, p :: Number, dt_f = 0.1, noise_factor = 2.0)
    g(u,p,t) = noise_factor*u*sin(u)
    prob = SDEProblem(prob_ode_linear!, g, u0, tspan, p)
    sol = solve(prob, saveat = dt_f)
    return sol.t, sol.u
end

# Solving the SDE problem using DifferentialEquations.jl with user given tspan
# u0, p are assigned (uniform) random values ∈ [-2.5,2.5]

SimpleLinearSDEModel(tspan :: Tuple{Number, Number}; dt_f = 0.1, noise_factor = 2.0) = SimpleLinearSDEModel(5*rand() - 2.5, tspan, 5*rand() - 2.5; dt_f = 0.1, noise_factor = 2.0*rand())


# Solving the ODE problem using DifferentialEquations.jl with user given u0, tspan
# p is assigned a (uniform) random values ∈ [-2.5,2.5]

SimpleLinearSDEModel(u0 :: Number, tspan :: Tuple{Number, Number}; dt_f = 0.1, noise_factor = 2.0) = SimpleLinearSDEModel(u0, tspan, typeof(u0)(5*rand() - 2.5), dt_f, noise_factor)

# Solving the ODE problem using DifferentialEquations.jl with user given tspan, p
# u0 is assigned a (uniform) random values ∈ [-2.5,2.5]

SimpleLinearSDEModel(tspan :: Tuple{Number, Number}, p :: Number; dt_f = 0.1, noise_factor = 2.0) = SimpleLinearSDEModel(typeof(p)(5*rand() - 2.5), tspan, p, dt_f, noise_factor)

# multiple dispatch copy for the original function

SimpleLinearSDEModel(u0 :: Number, tspan :: Tuple{Number, Number}, p :: Number; dt_f = 0.1, noise_factor = 2.0) = SimpleLinearSDEModel(u0, tspan, p, dt_f, noise_factor)


# ------------------------------------------------------------------------------


# example function call:
    # t, y = SimpleLinearODEModel(ode conditions; dt_f = ode time step)
    # t, y = SimpleLinearSDEModel(sde conditions; dt_f = ode time step, noise_factor = noise_factor)

testtime, testsol = SimpleLinearODEModel(1.0, (0.0, 2.0), 2.0; dt_f = 0.01)
scatter(testtime, testsol)
