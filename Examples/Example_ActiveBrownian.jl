using SoftSquishyMatter
using Random # useful if you want to shuffle arrays
cd(dirname(@__FILE__)) # set directory to current file location

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L_x = 100e-6 # width of simulation box
L_y = 100e-6 # height of simulation box

ptype = :activebrownian # particle label used to sorting and plotting
R = 1e-6 # radius of particle
γ_trans = 6 * pi * 8.9e-4 * R # friction coefficient of particle assuming Stokes flow
D_trans = 1.38e-23 * 300 / γ_trans # diffusivity of particle
φ = 0.4 # packing fraction

γv = γ_trans * 15e-6 # propulsion force of active Brownian particle
D_rot = 1.0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize simulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = Simulation() # create empty simulation

# description of simulation
simulation.descriptor = "Active Brownian particles with thermal noise"

# set simulation width and height
simulation.L_x = L_x
simulation.L_y = L_y

# initialize particles on a square lattice with packing fraction φ
for (x, y) in square_lattice_with_φ(φ = φ, L_x = L_x, L_y = L_y, radius = R)
    push!(simulation.particles, Particle(ptype = ptype, x = x, y = y, R = R, D_trans = D_trans, active_force = ActiveBrownian(γv = γv, D_rot = D_rot)))
end
pgroup_all = group_by_type(simulation.particles; ptype = :activebrownian)

# initialize cell list and lennard-Jones interaction
cell_list = CellList(particles = pgroup_all, L_x = L_x, L_y = L_y, cutoff = 2^(1.0 / 6.0) * 2 * R)
lj = LennardJones(particles = pgroup_all, cell_list = cell_list, ϵ = D_trans * γ_trans, multithreaded = true)
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
run_simulation(simulation; save_to = "out/Example_ActiveBrownian_data.out")

# OPTIONAL: for visualization and testing purposes
simulation = load_simulation(file = "out/Example_ActiveBrownian_data.out")
plot_history!(simulation.history; frame_size = (600, 600), xlim = [0.0, simulation.L_x], ylim = [0.0, simulation.L_y], colors = Dict(:activebrownian => "black"), folder = "frames/Example_ActiveBrownian/")
