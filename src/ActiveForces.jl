export get_active_force
export update_active_force

"""
    get_active_force(active_brownian)

Returns the components of the propulsion force `γv` for active Brownian motion.
This function is used to avoid some type stability issues.
"""
@inline function get_active_force(active_brownian::ActiveBrownian)
    return active_brownian.γv_x, active_brownian.γv_y
end

"""
    update_active_force!(active_brownian, particle; dt)

Advances the active Brownian force by a timestep `dt`.  Currently this
requires `particle` as an input, which shouldn't be necessary.  This will be
removed in future versions.
"""
function update_active_force(active_brownian::ActiveBrownian; particle_orientation::Float64, dt::Float64)
    if active_brownian.align
        active_brownian.θ = particle_orientation
    else
        active_brownian.θ += active_brownian.amplitude * sqrt(dt) * randn()
    end
    active_brownian.γv_x = active_brownian.γv * cos(active_brownian.θ)
    active_brownian.γv_y = active_brownian.γv * sin(active_brownian.θ)

    return particle_orientation
end

"""
    get_active_force(run_and_tumble)

Returns the components of the propulsion force `γv` for run-and-tumble
motion.  This function is used to avoid some type stability issues.
"""
@inline function get_active_force(run_and_tumble::RunAndTumble)
    return run_and_tumble.γv_x, run_and_tumble.γv_y
end

"""
    update_active_force!(active_brownian, particle; dt)

Advances the run-and-tumble force by a timestep `dt`.  Currently this
requires `particle` as an input, which shouldn't be necessary.  This will be
removed in future versions.
"""
function update_active_force(run_and_tumble::RunAndTumble; particle_orientation::Float64, dt::Float64)
    run_and_tumble.τ_run -= dt
    if run_and_tumble.τ_run < 0.0
        run_and_tumble.θ = 2 * pi * rand()
        run_and_tumble.τ_run = -log(rand()) / run_and_tumble.α

        if run_and_tumble.align
            particle_orientation = run_and_tumble.θ
        end

        run_and_tumble.γv_x = run_and_tumble.γv * cos(run_and_tumble.θ)
        run_and_tumble.γv_y = run_and_tumble.γv * sin(run_and_tumble.θ)
    else
        if run_and_tumble.align
             run_and_tumble.θ = particle_orientation

             run_and_tumble.γv_x = run_and_tumble.γv * cos(run_and_tumble.θ)
             run_and_tumble.γv_y = run_and_tumble.γv * sin(run_and_tumble.θ)
        end
    end
    # these cos and sin don't need to be computed every step.  Reorganize this function...


    return particle_orientation
end
