export update_particles!

"""
    update_particles!(brownian; period_x, period_y)

Advances the particles stored under `brownian.particles` by a timestep
`brownian.dt`.  Apply periodic boundaries with periodicity `period_x`,
`period_y` if `period_x > 0` or `period_y > 0`.
"""
function update_particles!(brownian::Brownian; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads brownian.multithreaded for particle in brownian.particles
        if !isnothing(particle.active_force)
            particle.θ = update_active_force(particle.active_force; particle_orientation = particle.θ, dt = brownian.dt)

            af_x, af_y = get_active_force(particle.active_force)
            particle.f_x += af_x
            particle.f_y += af_y
        end

        dtγ_trans = brownian.dt / particle.γ_trans
        @fastmath amplitude = sqrt(2 * particle.D_trans * brownian.dt)
        particle.x = wrap_position(particle.x + dtγ_trans * particle.f_x + amplitude * randn(); period = period_x)
        particle.y = wrap_position(particle.y + dtγ_trans * particle.f_y + amplitude * randn(); period = period_y)

        particle.f_x = 0.0
        particle.f_y = 0.0

        if brownian.rotations
            dtγ_rot = brownian.dt / particle.γ_rot
            @fastmath particle.θ += dtγ_rot * particle.t_θ + sqrt(2 * particle.D_rot * brownian.dt) * randn()

            particle.t_θ = 0.0
        end
    end
end
