module SoftSquishyMatter

println(args...) = println(stdout, args...)
function println(io::IO, args...)
    Base.println(io, args...)
    flush(io)
end

using Random
using Serialization
using DataStructures
using Plots

include("definitions.jl")
include("runsimulation.jl")
include("activeforces.jl")
include("particle.jl")
include("integrators.jl")
include("celllist.jl")
include("interactions.jl")
include("externalforces.jl")

include("analysis/visualize.jl")

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("SoftSquishyMatter Loaded!")
println("Michael Wang")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

end
