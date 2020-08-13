using SoftSquishyMatter # load simulation package
using Random # useful if you want to shuffle arrays
cd(dirname(@__FILE__)) # set directory to this file's folder

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize empty simulation with a description for future reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = Simulation()
simulation.descriptor = "Poly-dispersed colloids with Lennard-Jones interactions"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize hot and cold particles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R_small = 1e-6 # radius of particle
γ_trans_small = 6 * pi * 8.9e-4 * R_small # friction coefficient of particle assuming Stokes flow
D_trans_small = 1.38e-23 * 300 / γ_trans_small # diffusivity of particle

R_large = 5e-6 # radius of large particles
γ_trans_large = 6 * pi * 8.9e-4 * R_large
D_trans_large = 1.38e-23 * 300 / γ_trans_large

# initialize small and large particles on triangular lattice
# use the lattice dimensions of the small particles for the simulation dimensions (L_x, L_y)
positions_small, simulation.L_x, simulation.L_y = triangular_lattice(s = 3e-6, M_x = 30, M_y = 20)
positions_large, _, _ = triangular_lattice(s = 20e-6, M_x = 4, M_y = 3)

# remove small particles that overlap with the large ones
remove_overlaps!(positions_small; fixed_positions = positions_large, minimum_distance = 2^(1.0 / 6.0) * (R_small + R_large), period_x = simulation.L_x, period_y = simulation.L_y)

# create the particles with properties
for (x, y) in positions_small
    push!(simulation.particles, Particle(ptype = :small, x = x, y = y, R = R_small, γ_trans = γ_trans_small, D_trans = D_trans_small))
end
for (x, y) in positions_large
    push!(simulation.particles, Particle(ptype = :large, x = x, y = y, R = R_large, γ_trans = γ_trans_large, D_trans = D_trans_large))
end
# particle groups for small and large
# these will be useful for the interactions between small-small, small-large, and large-large
pgroup_small = group_by_type(simulation.particles; ptype = :small)
pgroup_large = group_by_type(simulation.particles; ptype = :large)
pgroup_all = simulation.particles

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize cell lists and Lennard-Jones interactions
# We here use different cell lists for the different particles.  This speeds up
# the simulation compared to if we used one cell list (see comment below).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cell_list_ss = CellList(particles = pgroup_small, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * 2 * R_small)
cell_list_sl = CellList(particles = pgroup_small, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * (R_small + R_large))
cell_list_ll = CellList(particles = pgroup_large, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * 2 * R_large)
lj_ss = LennardJones(particles = pgroup_small, cell_list = cell_list_ss, ϵ = 1.38e-23 * 300, multithreaded = true)
lj_sl = LennardJones(particles = pgroup_large, cell_list = cell_list_sl, ϵ = 1.38e-23 * 300, multithreaded = false, use_newton_3rd = true)
lj_ll = LennardJones(particles = pgroup_large, cell_list = cell_list_ll, ϵ = 1.38e-23 * 300, multithreaded = false)
append!(simulation.cell_lists, [cell_list_ss, cell_list_sl, cell_list_ll])
append!(simulation.interactions, [lj_ss, lj_sl, lj_ll])
# this works too, but is slower since small particles can use a smaller cutoff
#=cell_list = CellList(particles = pgroup_all, L_x = simulation.L_x, L_y = simulation.L_y, cutoff = 2^(1.0 / 6.0) * 2 * R_large)
lj = LennardJones(particles = pgroup_all, cell_list = cell_list, ϵ = 1.38e-23 * 300, multithreaded = true)
push!(simulation.cell_lists, cell_list)
push!(simulation.interactions, lj)=#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize Brownian integrator for overdamped dynamics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.dt = 0.001
# we apply the same integrator to all the particles
brownian = Brownian(particles = pgroup_all, dt = simulation.dt, multithreaded = true)
push!(simulation.integrators, brownian)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the number of steps and which particles to save periodically
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation.num_steps = trunc(Int64, 10 / simulation.dt)
simulation.save_interval = trunc(Int64, 1 / simulation.dt)
simulation.save_particles = pgroup_all

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run simulation and save data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_simulation!(simulation, save_as = "out/Example_PolyDispersed_data.out")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OPTIONAL: Load simulation data and animate the simulation.  Individual frames
# can be saved with plot_frames!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulation = load_simulation(file = "out/Example_PolyDispersed_data.out")
animate_frames!(simulation; colors = Dict(:small => "black", :large => "blue"), fps = 1, save_as = "frames/Example_PolyDispersed.gif")
