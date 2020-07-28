using SoftSquishyMatter # load simulation package
using Random # useful if you want to shuffle arrays
cd(dirname(@__FILE__)) # set directory to this file's folder

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize empty simulation with description for future reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = Simulation()
simulation.descriptor = "Run-and-tumble particles with thermal noise"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize particles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R = 1e-6 # radius of particle
γ_trans = 6 * pi * 8.9e-4 * R # friction coefficient of particle assuming Stokes flow
D_trans = 1.38e-23 * 300 / γ_trans # diffusivity of particle

# initialize positions of particles on a triangular lattice
# use lattice dimensions for the simulation dimensions (L_x, L_y)
positions, simulation.L_x, simulation.L_y = triangular_lattice(s = 3e-6, M_x = 30, M_y = 20)

# create particles with properties
for (x, y) in positions
    active_force = RunAndTumble(γv = γ_trans * 15e-6, α = 1.0)
    push!(simulation.particles, Particle(ptype = :runandtumble, x = x, y = y, R = R, γ_trans = γ_trans, D_trans = D_trans, active_force = active_force))
end
# particle groups are useful for interactions and integrators
# in this case, it doesn't really matter
pgroup_all = simulation.particles

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize cell list and Lennard-Jones interactions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cell_list = CellList(particles = pgroup_all, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * 2 * R)
lj = LennardJones(particles = pgroup_all, cell_list = cell_list, ϵ = 1.38e-23 * 300, multithreaded = true)
push!(simulation.cell_lists, cell_list)
push!(simulation.pair_interactions, lj)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize Brownian integrator for overdamped dynamics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.dt = 0.0001
brownian = Brownian(particles = pgroup_all, dt = simulation.dt, multithreaded = true)
push!(simulation.integrators, brownian)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the number of steps and which particles to save periodically
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.num_steps = trunc(Int64, 30 / simulation.dt)
simulation.save_interval = trunc(Int64, 1 / simulation.dt)
simulation.save_particles = pgroup_all

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run simulation and save data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_simulation(simulation; save_to = "out/Example_RunAndTumble_data.out")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OPTIONAL: Load simulation data and make each saved frame of simulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = load_simulation(file = "out/Example_RunAndTumble_data.out")
plot_history!(simulation.history; xlim = [0.0, simulation.L_x], ylim = [0.0, simulation.L_y], colors = Dict(:runandtumble => "black"), folder = "frames/Example_RunAndTumble/")
