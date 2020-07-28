using SoftSquishyMatter # load simulation package
using Random # useful if you want to shuffle arrays
cd(dirname(@__FILE__)) # set directory to this file's folder

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize empty simulation with a description for future reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = Simulation()
simulation.descriptor = "Run-and-tumble particles with funnels"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize active particles and boundary particles (walls and funnels)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R_swimmer = 1e-6 # radius of active particles
γ_trans_swimmer = 6 * pi * 8.9e-4 * R_swimmer # translational friction coefficient
D_trans_swimmer = 1.38e-23 * 300 / γ_trans_swimmer # translational diffusivity

R_boundary = 1e-6 # radius of boundary particles

# initialize positions of active particles on a triangular lattice
# use lattice dimensions for the simulation dimensions (L_x, L_y)
swimmer_positions, simulation.L_x, simulation.L_y = triangular_lattice(s = 4e-6, M_x = 50, M_y = 15)

# initialize boundary particles
boundary_positions = Array{Array{Float64, 1}, 1}()
# add two walls of particles on the leftmost and rightmost edges of simulation
append!(boundary_positions, flat_edge(from = [R_boundary, R_boundary], to = [R_boundary, simulation.L_y - R_boundary], spacing = R_boundary))
append!(boundary_positions, flat_edge(from = [simulation.L_x - R_boundary, R_boundary], to = [simulation.L_x - R_boundary, simulation.L_y - R_boundary], spacing = R_boundary))
# add lanes of wedges to form funnels
angle = pi / 2 # wedge angle
edge_length = 12e-6 # length of wedge arms
num_wedges_per_lane = 4 # number of wedges per lane
num_lanes = 3 # number of lanes of funnels
wedge_to_wedge = simulation.L_y / num_wedges_per_lane # distance between wedges
lane_to_lane = simulation.L_x / num_lanes # distance between lanes
for i = 0 : num_lanes - 1, j = 0 : num_wedges_per_lane - 1
    tip_position = [lane_to_lane / 2 + i * lane_to_lane - edge_length / (2 * sqrt(2)), wedge_to_wedge / 2 + j * wedge_to_wedge]
    append!(boundary_positions, wedge(angle = angle, edge_length = edge_length, spacing = R_boundary, tip = tip_position))
end

# remove active particles that overlap with boundary particles
remove_overlaps!(swimmer_positions; fixed_positions = boundary_positions, minimum_distance = 2^(1.0 / 6.0) * (R_swimmer + R_boundary), period_x = simulation.L_x, period_y = simulation.L_y)

# create the actual particles with properties
for (x, y) in swimmer_positions
    active_force = RunAndTumble(γv = γ_trans_swimmer * 15e-6, α = 1.0)
    push!(simulation.particles, Particle(ptype = :swimmer, x = x, y = y, R = R_swimmer, γ_trans = γ_trans_swimmer, D_trans = D_trans_swimmer, active_force = active_force))
end
for (x, y) in boundary_positions
    push!(simulation.particles, Particle(ptype = :boundary, x = x, y = y, R = R_boundary, D_trans = 0.0))
end
# particle groups are useful for interactions and integrators
pgroup_swimmer = group_by_type(simulation.particles; ptype = :swimmer)
pgroup_all = simulation.particles

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize cell list and Lennard-Jones interactions for swimmer-swimmer and swimmer-boundary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cell_list = CellList(particles = pgroup_all, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * (R_swimmer + R_boundary))
lj = LennardJones(particles = pgroup_swimmer, cell_list = cell_list, ϵ = 1.38e-23 * 300, multithreaded = true)
push!(simulation.cell_lists, cell_list)
push!(simulation.pair_interactions, lj)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize Brownian integrator for overdamped dynamics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.dt = 0.0001
# only swimmers will be integrated since boundary doesn't move
brownian = Brownian(particles = pgroup_swimmer, dt = simulation.dt, multithreaded = true)
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
run_simulation(simulation; save_to = "out/Example_RunAndTumbleFunnels_data.out")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OPTIONAL: Load simulation data and make each saved frame of simulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = load_simulation(file = "out/Example_RunAndTumbleFunnels_data.out")
plot_history!(simulation.history; xlim = [0.0, simulation.L_x], ylim = [0.0, simulation.L_y], colors = Dict(:boundary => "black", :swimmer => "green"), folder = "frames/Example_RunAndTumbleFunnels/")
