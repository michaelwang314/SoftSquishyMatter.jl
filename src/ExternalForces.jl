export compute_external_force!

#=
Constant force
=#
function compute_external_force!(constant_force::ConstantForce)
    for particle in constant_force.particles
        particle.f_x += constant_force.f_x
        particle.f_y += constant_force.f_y
    end
end

#=
Harmonic trap
=#
function compute_external_force!(harmonic_trap::HarmonicTrap)
    for particle in harmonic_trap.particles
        particle.f_x += -harmonic_trap.k_trap * (particle.x - harmonic_trap.x_center)
        particle.f_y += -harmonic_trap.k_trap * (particle.y - harmonic_trap.y_center)
    end
end
