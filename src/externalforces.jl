export compute_external_forces!

"""
    compute_external_force!(contant_force)

Applies a constant force to select particles stored under
`constant_force.particles`.
"""
function compute_external_forces!(constant_force::ConstantForce)
    @use_threads constant_force.multithreaded for particle in constant_force.particles
        particle.f_x += constant_force.f_x
        particle.f_y += constant_force.f_y
    end
end

"""
    compute_external_force!(harmonic_trap)

Applies a harmonic potential to select particles stored under
`harmonic_trap.particles`.
"""
function compute_external_forces!(harmonic_trap::HarmonicTrap)
    @use_threads harmonic_trap.multithreaded for particle in harmonic_trap.particles
        particle.f_x += -harmonic_trap.k_trap * (particle.x - harmonic_trap.x_center)
        particle.f_y += -harmonic_trap.k_trap * (particle.y - harmonic_trap.y_center)
    end
end
