using SoftSquishyMatter
using Random # useful if you want to shuffle arrays
cd(dirname(@__FILE__)) # set directory to current file location

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L_x = 100e-6 # width of simulation box
L_y = 100e-6 # height of simulation box

R = 1e-6 # radius of particle
γ_trans = 6 * pi * 8.9e-4 * R # friction coefficient of particle assuming Stokes flow
φ = 0.4 # packing fraction
fraction_hot = 0.4
ptype_hot = :hot # hot particle label
D_trans_hot = 1.38e-23 * 3000 / γ_trans # diffusivity of particle
ptype_cold = :cold # hot particle label
D_trans_cold = 1.38e-23 * 300 / γ_trans # diffusivity of particle

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize simulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = Simulation() # create empty simulation

# description of simulation
simulation.descriptor = "Mixture of Brownian particles at two different temperatures"

# set simulation width and height
simulation.L_x = L_x
simulation.L_y = L_y

# initialize particles on a square lattice with packing fraction φ
positions = square_lattice_with_φ(φ = φ, L_x = L_x, L_y = L_y, radius = R)
shuffle!(positions)
for (x, y) in positions[1 : trunc(Int64, length(positions) * fraction_hot)]
    push!(simulation.particles, Particle(ptype = ptype_hot, x = x, y = y, R = R, D_trans = D_trans_hot))
end
for (x, y) in positions[trunc(Int64, length(positions) * fraction_hot) + 1 : end]
    push!(simulation.particles, Particle(ptype = ptype_cold, x = x, y = y, R = R, D_trans = D_trans_cold))
end
pgroup_all = group_by_type(simulation.particles; ptype = [:cold, :hot])

# initialize cell list and lennard-Jones interaction
cell_list = CellList(particles = pgroup_all, L_x = L_x, L_y = L_y, cutoff = 2^(1.0 / 6.0) * 2 * R)
lj = LennardJones(particles = pgroup_all, cell_list = cell_list, ϵ = D_trans_hot * γ_trans, multithreaded = true)
push!(simulation.cell_lists, cell_list)
push!(simulation.pair_interactions, lj)

# initialzie Brownian integrator for overdamped dynamics
simulation.dt = 0.0001
brownian = Brownian(particles = pgroup_all, dt = simulation.dt, multithreaded = true)
push!(simulation.integrators, brownian)

# set total number of steps to run and which particles to store priodically
simulation.num_steps = trunc(Int64, 30 / simulation.dt)
simulation.save_interval = trunc(Int64, 0.1 / simulation.dt)
simulation.save_particles = pgroup_all

# run simulation and save data
run_simulation(simulation; save_to = "out/Example_TwoTemperature_data.out")

# OPTIONAL: for visualization and testing purposes
simulation = load_simulation(file = "out/Example_TwoTemperature_data.out")
plot_history!(simulation.history; frame_size = (600, 600), xlim = [0.0, simulation.L_x], ylim = [0.0, simulation.L_y], colors = Dict(ptype_cold => "blue", ptype_hot => "red"), folder = "frames/Example_TwoTemperature/")
