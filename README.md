# Parameter estimation models using Julia

This project was done as a part of a summer research internship conducted under the supervision of Dr. Yaroslav Molkov. The project was centered around creating various forward models to serve as datasets for parameter estimation models powered by Neural ODEs and genetic algorithm

# Folders 
## Forward_data

The forward data folder consists of various well known models like a simple linear model, the Fithugh-Nagumo model, LotkaVolterra systems, Lorenz oscillator, three body simulator and some more. The data produced can be catered to user-defined or randomized initial conditions with a way to opt for stochastic solutions. 

## Parameter Estimation

The parameter estimation folder consists of a genetic algorithm scheme illustrating the neural convergence of a randomly initialized LotkaVolterra solution to a given dataset produced using a pre-defined LotkaVolterra system. The file can be easily augmented to include other models from the forward_data folder. 

## Next updates (that will be made eventually, no realistic timeline here)

Easy goals:
1. The integration between the parameter estimation schemes and the forward data has to be made more accessible. 
2. More forward models need to be used. 

Hard Goals
1. Better optimization and personalization of training needs to be made.
2. Upgrade the network architecture and use more sophesticated methods than genetic learning
