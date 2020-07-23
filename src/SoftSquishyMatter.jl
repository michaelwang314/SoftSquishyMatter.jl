module SoftSquishyMatter

#=
Load files
=#
include("Definitions.jl")
include("RunSimulation.jl")
include("ActiveForces.jl")
include("Particle.jl")
include("Integrators.jl")
include("CellList.jl")
include("PairInteractions.jl")
include("ExternalForces.jl")

include("Analysis/Plotter.jl")

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Simulation Package Loaded")
println("Michael Wang, 07/21/2020")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

end
