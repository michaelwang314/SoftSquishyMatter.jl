using SoftSquishyMatter # load simulation package
using Random # useful if you want to shuffle arrays
cd(dirname(@__FILE__)) # set directory to this file's folder

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize empty simulation with a description for future reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = Simulation()
simulation.descriptor = "Mixture of Brownian particles at two different temperatures"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize particles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R = 1e-6 # radius of particle
fraction_hot = 0.4 # the fraction of particles that are hot
γ_trans_hot = 6 * pi * 8.9e-4 * R # friction coefficient
D_trans_hot = 1.38e-23 * 3000 / γ_trans_hot # diffusivity of hot particles
γ_trans_cold = 6 * pi * 8.9e-4 * R # friction coefficient
D_trans_cold = 1.38e-23 * 300 / γ_trans_cold # diffusivity of cold particles

# initialize positions of particles on a triangular lattice
# use lattice dimensions for the simulation dimensions (L_x, L_y)
positions, simulation.L_x, simulation.L_y = triangular_lattice(s = 3e-6, M_x = 30, M_y = 20)
# random positions so that hot and cold particles are uniformly mixed
shuffle!(positions)

# create actual particles with properties
for (x, y) in positions[1 : trunc(Int64, length(positions) * fraction_hot)]
    push!(simulation.particles, Particle(ptype = :hot, x = x, y = y, R = R, γ_trans = γ_trans_hot, D_trans = D_trans_hot))
end
for (x, y) in positions[trunc(Int64, length(positions) * fraction_hot) + 1 : end]
    push!(simulation.particles, Particle(ptype = :cold, x = x, y = y, R = R, γ_trans = γ_trans_cold, D_trans = D_trans_cold))
end
# particle groups are useful for interactions and integrators
# in this case, it doesn't really matter
pgroup_all = simulation.particles

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize cell list and Lennard-Jones interactions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cell_list = CellList(particles = pgroup_all, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * 2 * R)
lj = LennardJones(particles = pgroup_all, cell_list = cell_list, ϵ = D_trans_hot * γ_trans_hot, multithreaded = true)
push!(simulation.cell_lists, cell_list)
push!(simulation.interactions, lj)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize Brownian integrator for dynamics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.dt = 0.0001
brownian = Brownian(particles = pgroup_all, dt = simulation.dt, multithreaded = true)
push!(simulation.integrators, brownian)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the number of steps and which particles to save periodically
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.num_steps = trunc(Int64, 10 / simulation.dt)
simulation.save_interval = trunc(Int64, 1 / simulation.dt)
simulation.particles_to_save = pgroup_all

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run simulation and save data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_simulation!(simulation; save_as = "out/Example_TwoTemperature_data.out")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OPTIONAL: Load simulation data and animate the simulation.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = load_simulation(file = "out/Example_TwoTemperature_data.out")
visualize!(simulation; save_as = (:gif, "frames/Example_TwoTemperature.gif"), fps = 1, particle_colors = Dict(:cold => "blue", :hot => "red"))
