export update_particles!

#=
Apply periodicity
=#
@inline function wrap_position(position::Float64; period::Float64)
    if period > 0.0
        return mod(position, period)
    end
    return position
end

#=
Update particles with 'Brownian' integrator
=#
function update_particles!(integrator::Brownian; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads integrator.multithreaded for particle in integrator.particles
        if !isnothing(particle.active_force)
            update_active_force!(particle.active_force, particle; dt = integrator.dt)

            af_x, af_y = get_active_force(particle.active_force)
            particle.f_x += af_x
            particle.f_y += af_y
        end

        dtγ_trans = integrator.dt / particle.γ_trans
        amplitude = sqrt(2 * particle.D_trans * integrator.dt)
        particle.x = wrap_position(particle.x + dtγ_trans * particle.f_x + amplitude * randn(), period = period_x)
        particle.y = wrap_position(particle.y + dtγ_trans * particle.f_y + amplitude * randn(), period = period_y)

        if integrator.rotations
            dtγ_rot = integrator.dt / particle.γ_rot
            particle.θ += dtγ_rot * particle.t_θ + sqrt(2 * particle.D_rot * integrator.dt) * randn()
        end

        particle.f_x = 0.0
        particle.f_y = 0.0
        particle.t_θ = 0.0
    end
end
