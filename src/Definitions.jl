#=
################################################################################
Active/Self-propelling forces
    - Active Brownian
    - Run-and-tumble
################################################################################
=#
export AbstractActiveForce
export ActiveBrownian
export RunAndTumble

abstract type AbstractActiveForce end

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

#=
################################################################################
Particle
################################################################################
=#
export Particle

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

#=
################################################################################
CellList for finding neighbors
################################################################################
=#
export CellList

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

#=
################################################################################
Types of pair interactions
    - Lennard-Jones
################################################################################
=#
export AbstractPairInteraction
export LennardJones

abstract type AbstractPairInteraction end

struct LennardJones <: AbstractPairInteraction
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

struct DipoleDipole <: AbstractPairInteraction
    #to be added
end

struct HarmonicBond <: AbstractPairInteraction
    #to be added
end

#=
################################################################################
External forces
################################################################################
=#
export AbstractExternalForce

abstract type AbstractExternalForce end

struct ConstantForce <: AbstractExternalForce
    f_x::Float64
    f_y::Float64

    particles::Array{Particle, 1}

    function ConstantForce(particles::Array{Particle, 1}; f_x::Float64, f_y::Float64)
        new(f_x, f_y, particles)
    end
end

struct HarmonicTrap <: AbstractExternalForce
    x_center::Float64
    y_center::Float64
    k_trap::Float64

    particles::Array{Particle, 1}

    function HarmonicTrap(particles::Array{Particle, 1}; x_center::Float64, y_center::Float64, k_trap::Float64)
        new(x_center, y_center, k_trap, particles)
    end
end

#=
################################################################################
Integrators
    - Brownian
################################################################################
=#
export AbstractIntegrator
export Brownian

abstract type AbstractIntegrator end

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

#=
################################################################################
For storing the simulation
################################################################################
=#
export Simulation

mutable struct Simulation
    descriptor::String

    L_x::Float64
    L_y::Float64
    periodic_in_x::Bool
    periodic_in_y::Bool

    particles::Array{Particle, 1}

    cell_lists::Array{CellList, 1}
    pair_interactions::Array{AbstractPairInteraction, 1}
    external_forces::Array{AbstractExternalForce, 1}

    dt::Float64
    integrators::Array{AbstractIntegrator, 1}

    num_steps::Int64
    save_interval::Int64
    save_particles::Array{Particle, 1}
    overwrite::Bool
    history::Array{Array{Particle, 1}}

    function Simulation()
        new("No description given...",
            0.0, 0.0, true, true,
            Array{Particle, 1}(),
            Array{CellList, 1}(), Array{AbstractPairInteraction, 1}(), Array{AbstractExternalForce, 1}(),
            0.0, Array{AbstractIntegrator, 1}(),
            0, 0, Array{Particle, 1}(), true, Array{Array{Particle, 1}, 1}())
    end
end
