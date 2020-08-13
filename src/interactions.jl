export compute_interactions!

const LJ_CONST = 1.2599210498948732 # 2^(1.0 / 3.0)

"""
    lj_coef(ϵ, σ², Δr²)

Returns ϵ / σ² * [48 * (σ² / Δr²)³ - 24] * (σ² / Δr²)⁴.  Useful for Lennard-
Jones and Shifted Lennard Jones.
"""
@inline function lj_coef(ϵ::Float64, σ²::Float64, Δr²::Float64)
    inv² = σ² / Δr²
    inv⁶ = inv² * inv² * inv²
    inv⁸ = inv⁶ * inv²
    return ϵ / σ² * (48 * inv⁶ - 24) * inv⁸
end

"""
    compute_interations!(lj; period_x, period_y)

Computes a Lennard-Jones force for particles stored under `lj.particles` given
neighbors stored in the cell list `lj.cell_list`.  If `period_x > 0` and
`period_y > 0`, then apply periodic boundary conditions to the pair interaction.
"""
function compute_interactions!(lj::LennardJones; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads lj.multithreaded for particle in lj.particles
        x, y = particle.x, particle.y
        i = trunc(Int64, x / lj.cell_list.cell_spacing_x)
        j = trunc(Int64, y / lj.cell_list.cell_spacing_y)
        f_x, f_y = 0.0, 0.0
        @inbounds for di = -1 : 1, dj = -1 : 1
            idi = mod(i + di, lj.cell_list.num_cells_x) + 1
            jdj = mod(j + dj, lj.cell_list.num_cells_y) + 1
            pid = lj.cell_list.start_pid[idi, jdj]
            while pid > 0
                neighbor = lj.cell_list.particles[pid]
                Δx = wrap_displacement(x - neighbor.x; period = period_x)
                Δy = wrap_displacement(y - neighbor.y; period = period_y)
                Δr² = Δx^2 + Δy^2

                σ² = (lj.σ < 0.0 ? particle.R + neighbor.R : lj.σ)^2
                cutoff² = (lj.cutoff < 0.0 ? LJ_CONST * σ² : lj.cutoff^2)
                if 0.0 < Δr² < cutoff²
                    coef = lj_coef(lj.ϵ, σ², Δr²)
                    f_x += (_f_x = coef * Δx)
                    f_y += (_f_y = coef * Δy)
                    if lj.use_newton_3rd
                        neighbor.f_x -= _f_x
                        neighbor.f_y -= _f_y
                    end
                end
                pid = lj.cell_list.next_pid[pid]
            end
        end
        particle.f_x += f_x
        particle.f_y += f_y
    end
end

"""
    compute_interations!(hbond; period_x, period_y)

Compute a harmonic force between pairs of particles stored under
`hb.particle_pairs`.  Apply periodic boundary conditions if `period_x > 0` or
`period_y > 0`.
"""
function compute_interactions!(hbond::HarmonicBond; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads hbond.multithreaded for (particle_1, particle_2) in hbond.particle_pairs
        Δx = wrap_displacement(particle_1.x - particle_2.x; period = period_x)
        Δy = wrap_displacement(particle_1.y - particle_2.y; period = period_y)
        Δr = sqrt(Δx^2 + Δy^2)

        # this check isn't really necessary unless user starts the particles in the same place
        if Δr > 0.0
            coef = -hbond.k_bond * (1 - hbond.l_rest / Δr)
            f_x = coef * Δx
            f_y = coef * Δy

            particle_1.f_x += f_x
            particle_1.f_y += f_y
            particle_2.f_x -= f_x
            particle_2.f_y -= f_y
        end
    end
end

"""
    compute_interations!(hθ; period_x, period_y)
"""
function compute_interactions!(hθ::HarmonicAngle; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads hθ.multithreaded for (particle_left, particle_mid, particle_right) in hθ.particle_triplets
        0
    end
end

"""
    compute_interations!(hcosθ; period_x, period_y)
"""
function compute_interactions!(hcosθ::HarmonicCosineAngle; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads hcosθ.multithreaded for (particle_left, particle_mid, particle_right) in hcosθ.particle_triplets
        0
    end
end
