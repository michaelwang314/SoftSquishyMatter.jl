export get_active_force
export update_active_force!

"""
    get_active_force(active_brownian)

Returns the x, y components of the propulsion force γv.  This function
"""
@inline function get_active_force(active_brownian::ActiveBrownian)
    return active_brownian.γv_x, active_brownian.γv_y
end

function update_active_force!(active_brownian::ActiveBrownian, particle::Particle; dt::Float64)
    if active_brownian.align
        active_brownian.θ = particle.θ
    else
        active_brownian.θ += active_brownian.amplitude * sqrt(dt) * randn()
    end
    active_brownian.γv_x = active_brownian.γv * cos(active_brownian.θ)
    active_brownian.γv_y = active_brownian.γv * sin(active_brownian.θ)
end

#=
Run-and-tumble
=#
@inline function get_active_force(run_and_tumble::RunAndTumble)
    return run_and_tumble.γv_x, run_and_tumble.γv_y
end

function update_active_force!(run_and_tumble::RunAndTumble, particle::Particle; dt::Float64)
    run_and_tumble.τ_run -= dt
    if run_and_tumble.τ_run < 0.0
        run_and_tumble.θ = 2 * pi * rand()
        run_and_tumble.τ_run = -log(rand()) / run_and_tumble.α

        if run_and_tumble.align
            particle.θ = run_and_tumble.θ
        end
    else
        if run_and_tumble.align
             run_and_tumble.θ = particle.θ
        end
    end
    run_and_tumble.γv_x = run_and_tumble.γv * cos(run_and_tumble.θ)
    run_and_tumble.γv_y = run_and_tumble.γv * sin(run_and_tumble.θ)
end
