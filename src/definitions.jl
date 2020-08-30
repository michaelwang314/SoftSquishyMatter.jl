#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Active forces
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractActiveForce
export ActiveBrownian
export RunAndTumble

abstract type AbstractActiveForce end

"""
Stores properties of active Brownian force

    ActiveBrownian(; γv, D_rot, align)

Initializes an active Brownian force with propulsion `γv`, rotational diffusion
`D_rot`.  If `align == true`, the active Brownian force will align with the
particle's orientation and allow us to introduce orientational interactions.
"""
mutable struct ActiveBrownian <: AbstractActiveForce
    γv_x::Float64
    γv_y::Float64
    θ::Float64

    γv::Float64
    amplitude::Float64

    align::Bool

    function ActiveBrownian(; γv::Float64, D_rot::Float64, align::Bool = false)
        θ = 2 * pi * rand()
        amplitude = sqrt(2 * D_rot)
        new(γv * cos(θ), γv * sin(θ), θ, γv, amplitude, align)
    end
end

"""
Stores properties of run-and-tumble force

    RunAndTumble(; γv, α, align)

Initializes a run-and-tumble force with propulsion `γv` and tumble rate `α`. If
`align == true`, the active Brownian force will align with the particle's
orientation and allow us to introduce orientational interactions.
"""
mutable struct RunAndTumble <: AbstractActiveForce
    γv_x::Float64
    γv_y::Float64
    θ::Float64
    τ_run::Float64

    γv::Float64
    α::Float64

    align::Bool

    function RunAndTumble(; γv::Float64, α::Float64, align::Bool = false)
        θ = 2 * pi * rand()
        τ_run = -log(rand()) / α

        new(γv * cos(θ), γv * sin(θ), θ, τ_run, γv, α, align)
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Particle
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export Particle

"""
Stores properties of a particle

    Particle(; ptype, x, y, θ, R, γ_trans, γ_rot, D_trans, D_rot, active_force)

# Arguments
- `ptype::Symbol = :particle`: particle type
- `x::Float64`: x position
- `y::Float64`: y position
- `θ::Float64 = 2 * pi * rand()`: in-plane orientation
- `R::Float64`: radius
- `γ_trans::Float64 = -1.0`: translational friction.  If `γ_trans < 0`, use `γ_trans = 6 * pi * 8.9e-4 * R`.
- `γ_rot::Float64 = -1.0`: rotational friction. If `γ_rot < 0`, use `γ_rot = 8 * pi * 8.9e-4 * R^3`.
- `D_trans::Float64`: translational diffusion
- `D_rot::Float64`: rotational diffusion
- `active_force::Union{AbstractActiveForce, Nothing} = nothing`: active force.
"""
mutable struct Particle
    ptype::Symbol # particle type

    x::Float64 # x position
    y::Float64 # y position
    θ::Float64 # orientation

    f_x::Float64 # x force
    f_y::Float64 # y force
    t_θ::Float64 # torque in z

    R::Float64 # particle radius
    γ_trans::Float64 # translational friction
    γ_rot::Float64 # rotational friction
    D_trans::Float64 # translational diffusivity
    D_rot::Float64 # rotational diffusivity

    active_force::Union{AbstractActiveForce, Nothing} # active force/self-propulsion

    function Particle(; ptype::Symbol = :particle,
                        x::Float64, y::Float64, θ::Float64 = 2 * pi * rand(),
                        R::Float64, γ_trans::Float64 = -1.0, γ_rot::Float64 = -1.0, D_trans::Float64, D_rot::Float64 = 0.0,
                        active_force::Union{AbstractActiveForce, Nothing} = nothing)
        if γ_trans < 0.0
            γ_trans = 6 * pi * 8.9e-4 * R
        end
        if γ_rot < 0.0
            γ_rot = 8 * pi * 8.9e-4 * R^3
        end

        new(ptype, x, y, θ, 0.0, 0.0, 0.0, R, γ_trans, γ_rot, D_trans, D_rot, active_force)
    end
end

"""
    Molecule(particles, pairs, triplets)

...
"""
struct Molecule
    particles::Array{Particle, 1}
    pairs::Array{NTuple{2, Particle}, 1}
    triplets::Array{NTuple{3, Particle}, 1}

    function Molecule(; particles::Array{Particle, 1} = Array{Particle, 1}(), 
                        pairs::Array{NTuple{2, Particle}, 1} = Array{NTuple{2, Particle}, 1}(), 
                        triplets::Array{NTuple{3, Particle}, 1} = Array{NTuple{3, Particle}, 1}())
        new(particles, pairs, triplets)
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cell list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export CellList

"""
Stores properties of a cell list

    CellList(; particles, L_x, L_y, cutoff)

Creates a cell list for `particles` in region/simulation with dimensions
`L_x`, `L_y`.  `cutoff` gives the approximate size of the cells, the actual size
of which is chosen so that an integer number of cells fit in each dimension.
"""
struct CellList
    start_pid::Array{Int64, 2}
    next_pid::Array{Int64, 1}

    particles::Array{Particle, 1}

    num_cells_x::Int64
    num_cells_y::Int64
    cell_spacing_x::Float64
    cell_spacing_y::Float64

    function CellList(; particles::Array{Particle, 1}, L_x::Float64, L_y::Float64, cutoff::Float64)
        num_cells_x = trunc(Int64, L_x / cutoff)
        num_cells_y = trunc(Int64, L_y / cutoff)
        cell_spacing_x = L_x / num_cells_x
        cell_spacing_y = L_y / num_cells_y

        start_pid = -ones(Int64, num_cells_x, num_cells_y)
        next_pid = -ones(Int64, length(particles))

        for (n, particle) in enumerate(particles)
            i = trunc(Int64, particle.x / cell_spacing_x) + 1
            j = trunc(Int64, particle.y / cell_spacing_y) + 1

            if start_pid[i, j] > 0
                next_pid[n] = start_pid[i, j]
            end
            start_pid[i, j] = n
        end

        new(start_pid, next_pid, particles, num_cells_x, num_cells_y, cell_spacing_x, cell_spacing_y)
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interactions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractInteraction
export LennardJones
export HarmonicBond
export HarmonicCosineAngle

abstract type AbstractInteraction end

"""
Stores properties of a Lennard-Jones interaction

    LennardJones(; particles, cell_list, ϵ, σ, cutoff, multithreaded, use_newton_3rd)
...
# Arguments
- `particles::Array{Particle, 1}`: particles to compute Lennard-Jones for
- `cell_list::CellList`: cell list for neighbors
- `ϵ::Float64`: energy scale of Lennard-Jones potential
- `σ::Float64 = -1.0`: length scale of Lennard-Jones potential.  If `σ < 0`, use diameters of interacting particles.
- `cutoff::Float64 = -1.0`: distance to truncate interaction.  If `cutoff < 0`, use `cutoff = 2^(1.0 / 6.0) * σ`, which corresponds to WCA interaction.
- `multithreaded::Bool = false`: if `true`, split particles between threads.
- `use_newton_3rd::Bool = false`: if `true`, use Newton's 3rd law to also compute force on neighbors.  Be CAREFUL with overcounting and race conditions with threads.
...
"""
struct LennardJones <: AbstractInteraction
    ϵ::Float64
    σ::Float64
    cutoff::Float64

    particles::Array{Particle, 1}
    cell_list::CellList

    multithreaded::Bool
    use_newton_3rd::Bool

    function LennardJones(; particles::Array{Particle, 1}, cell_list::CellList,
                            ϵ::Float64, σ::Float64 = -1.0, cutoff::Float64 = -1.0,
                            multithreaded::Bool = false, use_newton_3rd::Bool = false)
        new(ϵ, σ, cutoff, particles, cell_list, multithreaded, use_newton_3rd)
    end
end

"""
Stores properties of a harmonic bond between pair of particles

    HarmonicBond(; pairs, k_bond, l_rest, multithreaded)

Initialize a harmonic bond between `pairs` of particles with stiffness `k_bond`
and rest length `l_rest`.  For now `multithreaded` should only really be used 
for dimers to avoid race conditions.  This will be changed with the
introduction of Molecules.
"""
struct HarmonicBond <: AbstractInteraction
    k_bond::Float64
    l_rest::Float64

    particle_pairs::Array{NTuple{2, Particle}, 1}

    multithreaded::Bool

    function HarmonicBond(; pairs::Array{NTuple{2, Particle}, 1}, k_bond::Float64, l_rest::Float64 = 0.0, multithreaded::Bool = false)
        new(k_bond, l_rest, pairs, multithreaded)
    end
end

"""
    etc

etc
"""
struct HarmonicAngle <: AbstractInteraction
    k_θ::Float64
    θ_rest::Float64

    particle_triplets::Array{NTuple{3, Particle}, 1}

    multithreaded::Bool

    function HarmonicAngle(; triplets::Array{NTuple{3, Particle}, 1}, k_θ::Float64, θ_rest::Float64, multithreaded::Bool = false)
        new(k_θ, θ_rest, triplets, multithreaded)
    end
end

"""
    ...
    
...
"""
struct CosineAngle <: AbstractInteraction

end

"""
Stores properties of bending force for triplets of particles

    HarmonicCosineAngle(triplets, k_cosθ, cosθ_rest, multithreaded)

Initialize a bending force for `triplets` of particles with 
with bending stiffness `k_cosθ` with rest cosine angle 
`cosθ_rest`.  For now `multithreaded` should only really be used 
for trimers to avoid race conditions.  This will be changed with 
the introduction of Molecules.
"""
struct HarmonicCosineAngle <: AbstractInteraction
    k_cosθ::Float64
    cosθ_rest::Float64

    particle_triplets::Array{NTuple{3, Particle}, 1}

    multithreaded::Bool

    function HarmonicCosineAngle(; triplets::Array{NTuple{3, Particle}, 1}, k_cosθ::Float64, cosθ_rest::Float64, multithreaded::Bool = false)
        new(k_cosθ, cosθ_rest, triplets, multithreaded)
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# External forces
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractExternalForce

abstract type AbstractExternalForce end

"""
Stores properties of constant force

    ConstantForce(particles; f_x, f_y)

Initialize a constant force `f_x`, `f_y` for 'particles'.
"""
struct ConstantForce <: AbstractExternalForce
    f_x::Float64
    f_y::Float64

    particles::Array{Particle, 1}

    multithreaded::Bool

    function ConstantForce(particles::Array{Particle, 1}; f_x::Float64, f_y::Float64, multithreaded::Bool = false)
        new(f_x, f_y, particles, multithreaded)
    end
end

"""
Stores properties of harmonic trap

    HarmonicTrap(particles; x_center, y_center, k_trap)

Initialize a harmonic trap with stiffness `k_trap` centered at `x_center`,
`y_center` for `particles`.
"""
struct HarmonicTrap <: AbstractExternalForce
    x_center::Float64
    y_center::Float64
    k_trap::Float64

    particles::Array{Particle, 1}

    multithreaded::Bool

    function HarmonicTrap(particles::Array{Particle, 1}; x_center::Float64, y_center::Float64, k_trap::Float64, multithreaded::Bool = false)
        new(x_center, y_center, k_trap, particles, multithreaded)
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Integrators
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractIntegrator
export Brownian

abstract type AbstractIntegrator end

"""
Stores properties of a Brownian integrator

    Brownian(; particles, dt, rotations, multithreaded)

Initialize a Brownian integrator for `particles` with timestep `dt`.  If
`rotations == true`, integrate the orientational degree of freedom.  If
`multithreaded == true`, split particles between threads.
"""
struct Brownian <: AbstractIntegrator
    dt::Float64
    particles::Array{Particle, 1}

    rotations::Bool

    multithreaded::Bool

    function Brownian(; particles::Array{Particle, 1}, dt::Float64, rotations::Bool = false, multithreaded::Bool = false)
        new(dt, particles, rotations, multithreaded)
    end
end

struct Langevin <: AbstractIntegrator
    # to be added
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export Simulation

"""
Stores all information needed to run simulation

    Simulation(...)
...
# Parameters
- `descriptor::String`: a description of the simulation for easy referencing
- `L_x::Float64`: width of simulation region
- `L_y::Float64`: height of simulation region
- `periodic_in_x::Bool = true`: if `true`, apply periodicity along x
- `periodic_in_y::Bool = true`: if `true`, apply periodicity along y
- `particles::Array{Particle, 1}`: all particles in simulation
- `cell_lists::Array{CellList, 1}`: all cell lists used
- `interactions::Array{AbstractInteraction, 1}`: all interactions used
- `external_forces::Array{AbstractExternalForce, 1}`: all external forces used
- `dt::Float64`: timestep
- `integrators::Array{AbstractIntegrator}`: all integrators used
- `num_steps::Int64`: total steps to run
- `save_interval::Int64`: store data periodically
- `save_particles::Array{Particle, 1}`: particles to store
- `history::Array{Array{Particle, 1}}`: history of particles
...
"""
mutable struct Simulation
    descriptor::String

    L_x::Float64
    L_y::Float64
    periodic_in_x::Bool
    periodic_in_y::Bool

    particles::Array{Particle, 1}
    bonds::Array{NTuple{2, Particle}, 1}
    angles::Array{NTuple{3, Particle}, 1}

    cell_lists::Array{CellList, 1}
    interactions::Array{AbstractInteraction, 1}
    external_forces::Array{AbstractExternalForce, 1}

    dt::Float64
    integrators::Array{AbstractIntegrator, 1}
    num_steps::Int64

    save_interval::Int64
    things_to_save::NamedTuple{(:particles, :bonds), Tuple{Array{Particle, 1}, Array{NTuple{2, Particle}, 1}}}
    history::Array{NamedTuple{(:particles, :bonds), Tuple{Array{Particle, 1}, Array{NTuple{2, Particle}, 1}}}, 1}

    function Simulation()
        new("No description given...",
            0.0, 0.0, true, true,
            [], [], [],
            [], [], [],
            0.0, [], 0, 
            0, (particles = [], bonds = []), [])
    end
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Useful functions and macros
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export use_threads
export wrap_displacement
export wrap_positions
export print_message

"""
    @use_threads multithreaded ...

Applies `Threads.@threads` if `multithreaded == true`.  This is mostly used to
shorten code.
"""
macro use_threads(multithreaded::Union{Expr, Symbol}, expr::Expr)
    esc(quote
        if $multithreaded
            Threads.@threads $expr
        else
            $expr
        end
    end)
end

"""
    wrap_displacement(displacement; period)

Returns a new displacement after applying periodic boundary conditions.  The
periodicity is given by `period`.  If period < 0, then no periodicity is
applied.
"""
@inline function wrap_displacement(displacement::Float64; period::Float64)
    if period > 0.0 && abs(displacement) > period / 2
        return displacement - sign(displacement) * period
    end
    return displacement
end

"""
    wrap_position(position; period)

Returns a new position after applying periodic boundary conditions.  The
periodicity is given by `period`.  If `period < 0`, then no periodic boundary
condition is applied.
"""
@inline function wrap_position(position::Float64; period::Float64)
    if period > 0.0
        return mod(position, period)
    end
    return position
end

"""

"""
function print_message(message::String)
    message = string("++++++ ", message, " ++++++")
    N_space = max(div(62 - length(message), 2), 0)
    for _ = 1 : N_space
        message = string(" ", message)
    end

    println("")
    println(message)
    println("")
end
