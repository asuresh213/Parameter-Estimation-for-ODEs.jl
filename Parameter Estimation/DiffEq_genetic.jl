using DifferentialEquations
using DiffEqSensitivity
using Plots
using Flux

# Learning to fit a diff-eq solution using genetic algorithm + learning
# each generation has 5 members, whose genes are the parameters.
# in each generation, the members are given 500 epochs to match the data
# at the end of the first generation, the member with the lowest loss is taken
# its genes are replicated by copying them into the next generation and adding a small bias
# the offsprings now equiped with the superior genes converge to the solution faster than their previous generation


#---------------------------- defining the model -------------------------------

lv! = @ode_def begin
    dx = a*x - b*x*y
    dy = -c*y + d*x*y
end a b c d

u0 = [1.0, 1.0]
p0 = (2.5, 3.0, 2.0, 1.0) #adjust as per your liking
tspan = (0.0, 10.0)
tsteps = 0.0:0.1:10.0
#-------------------------------------------------------------------------------

#------------------- making the data that the algorithm will fit ---------------

function make_data(f, u0, p0, tspan)
    prob0 = ODEProblem(lv!, u0, tspan, p0)
    sol0 = solve(prob0, Tsit5(), saveat = 0.1)

    #creating discrete data!
    prey = sol0[1, :]
    pred = sol0[2, :]

    # plotting the data
    gr()
    scatter(tsteps, prey, label = "prey")
    scatter!(tsteps, pred, label = "predator")
    return prob0, sol0, prey, pred
end

prob0, sol0, prey, pred = make_data(lv!, u0, p0, tspan)


# Note: The above step can be replaced with the data from any of the ode systems
# from the forward models folder. However, one needs to (for the time being)
# re-write the ode system here and declare it to be prob0 for everything to work
# smoothly.

# Creating a nice wrapper for this shouldn't be to hard, and it is coming up in
# my to-do list.


# uncomment the following function to introduce noise into the data.
# Note: Unfortunately, introducing noise throws off the exit condition for the
# genetic model. So refactoring required!
# maybe try an SDE problem instead, with noise added in?

#=function noise!(sol0)
    noise1 = 0.5*rand(length(sol0))
    noise2 = 0.5*rand(length(sol0))
    sol0[1, :] += noise1
    sol0[2, :] += noise2
end
noise!(sol0)=#

#------------------------- Flux training parameters ----------------------------

# defining the loss function
function loss()
    prediction = concrete_solve(prob0, Tsit5(), u0, p, saveat = 0.1, sensealg = TrackerAdjoint())
    return sum(abs2, prediction - sol0)
end

# callback function to visualize training
iter = 0 # dummy variable for cb
cb = function()
    global iter += 1
    if (iter%10 == 1)
        scatter(tsteps, prey, label = "prey data")
        scatter!(tsteps, pred, label = "pred data")
        remade_sol = solve(remake(prob0, p=p), Tsit5(), saveat = 0.1)
        remade_prey = remade_sol[1, :]
        remade_pred = remade_sol[2, :]
        plot!(tsteps, remade_prey, ylim = (0, 6), labels = "prey model", color = 1)
        display(plot!(tsteps, remade_pred, ylim = (0, 6), labels = "pred model", color = 2))
    end
end

opt = ADAM(0.1) # optimizer for the descent
ldata = Iterators.repeated((),500) # defining the number of epochs


#-------------------- Genetic Algorithm helper functions -----------------------



# finding the best member of the population

function find_best(population, max_score)
    for ind in population
        if(ind.score == max_score)
            best = ind
            return best
        end
    end
end


# mutating a randomly picked induvidual

function mutation!(rand_pick)
    rand_num = rand(1:90)
    if(rand_num <= 30)
        rand_pick.gene += 3*rand(4) - rand(4)
        rand_pick.gene = abs.(rand_pick.gene)
    elseif(rand_num >= 30 && rand_num < 60)
        rand_pick.gene += 2*rand(4) .- 0.5*rand(4)
        rand_pick.gene = abs.(rand_pick.gene)
    elseif(rand_num >= 60 && rand_num < 90)
        rand_pick.gene += rand(4) .- 0.15*rand(4)
        rand_pick.gene = abs.(rand_pick.gene)
    end
end



#-------------------- Setting up the genetic algorithm -------------------------

max_population = 5 # maximum population per generation

# create a mutable struct for member properties
mutable struct member
    id
    gene
    score
end


# give your best guess for parameters for a good start point.
    # Default: [0.0, 0.0, 0.0, 0.0] for bare-bones genetic algorithm.
# Note: Default options might have very slow convergence for problems
# with bigger params

guess_param = [0.0, 0.0, 0.0, 0.0]


# create the first population

if(guess_param == [0.0, 0.0, 0.0, 0.0])
    population = [member(i, 1.2*rand(4) .+ 0.01, 0) for i in 1:max_population]
else
    population = [member(i, guess_param + 0.01*rand(4), 0) for i in 1:max_population]
end

# training each member of the population

fin = false # to keep track of convergence
gen = 1 # initial generation number

p = [0, 0, 0, 0] # initial params vector (will be changed while training)



while(fin == false)
    # ---- training phase ----
    @info "Training generation $gen"
    for ind in population
        global p = ind.gene # accessing the gene as the parameter
        cb() # visualize the initial set up
        @info "--Training induvidual $(ind.id)"
        ps = params(p) # params for the flux model
        Flux.train!(loss, ps, ldata, opt, cb = Flux.throttle(cb, 0.2)) # training
        ind.score = 1/loss() # scoring each induvidual based on their loss

        # exit condition (if loss < 0.01)
        if(σ(log(ind.score)) > 0.995)
            global fin = true
            println("The parameters of the given model is $(ind.gene)")
            @info "thanks to member $(ind.id) from generation $gen"
            break
        end

    end

    if(fin == false)
        # ---- replicating dominant genes phase ----
        max_score = maximum([ind.score for ind in population]) # finding maximum score from generation
        # finding the best induvidual based on score
        best = find_best(population, max_score) #best member of generation
        dom_gene = best.gene #copying the best genes into dominant gene
        # alter population to produce offsprings with genes = (dominant gene ± U[-0.05, 0.05]) and base score = 0.
        global population = [member(i, abs.(dom_gene + (0.1*rand(4) .- 0.05)), 0) for i in 1:max_population]

        # ---- mutation phase ----
        # pick a random member of the population to mutate
        rand_pick = rand(population)
        mutate = rand()
        if(mutate > 0.5)
            mutation!(rand_pick)
        end

        # ---- passing to next generation ----
        # increase generation count
        global gen +=1
        # if it takes too long, exit!
        if(gen >= 10)
            println("Convergence took too long.")
            global fin = true
        end
    end
end
