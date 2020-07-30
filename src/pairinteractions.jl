export compute_pair_interaction!

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
    compute_pair_interaction!(lj)

Computes a Lennard-Jones force for particles stored under `lj.particles` given
neighbors stored in the cell list `lj.cell_list`.  If `period_x > 0` and
`period_y > 0`, then apply periodic boundary conditions to the pair interaction.
"""
function compute_pair_interaction!(lj::LennardJones; period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    @use_threads lj.multithreaded for particle in lj.particles
        x = particle.x
        y = particle.y
        i_box = trunc(Int64, x / lj.cell_list.cell_spacing_x)
        j_box = trunc(Int64, y / lj.cell_list.cell_spacing_y)
        for di = -1 : 1, dj = -1 : 1
            i = mod(i_box + di, lj.cell_list.num_cells_x) + 1
            j = mod(j_box + dj, lj.cell_list.num_cells_y) + 1

            pid = lj.cell_list.start_pid[i, j]
            while pid > 0
                neighbor = lj.cell_list.particles[pid]
                Δx = wrap_displacement(x - neighbor.x; period = period_x)
                Δy = wrap_displacement(y - neighbor.y; period = period_y)

                Δr² = Δx^2 + Δy^2
                σ² = (lj.σ < 0.0 ? particle.R + neighbor.R : lj.σ)^2
                cutoff² = (lj.cutoff < 0.0 ? 2^(1.0 / 3.0) * σ² : lj.cutoff^2)
                if 0.0 < Δr² < cutoff²
                    inv² = σ² / Δr²
                    inv⁶ = inv² * inv² * inv²
                    inv⁸ = inv⁶ * inv²
                    coef = 24 * lj.ϵ / σ² * (2 * inv⁶ - 1) * inv⁸

                    f_x = coef * Δx
                    f_y = coef * Δy
                    particle.f_x += f_x
                    particle.f_y += f_y

                    if lj.use_newton_3rd
                        neighbor.f_x -= f_x
                        neighbor.f_y -= f_y
                    end
                end
                pid = lj.cell_list.next_pid[pid]
            end
        end
    end
end

"""
    compute_pair_interaction!(hbond; period_x, period_y)

Compute a harmonic force between pairs of particles stored under
`hb.particle_pairs`.  Apply periodic boundary conditions if `period_x > 0` or
`period_y > 0`
"""
function compute_pair_interaction!(hb::HarmonicBond; period_x = -1.0, period_y = -1.0)
    @use_threads hb.multithreaded for (particle1, particle2) in hb.particle_pairs
        x1, y1 = particle1.x, particle1.y
        x2, y2 = particle2.x, particle2.y

        Δx = wrap_displacement(x1 - x2; period = period_x)
        Δy = wrap_displacement(y1 - y2; period = period_y)
        Δr = sqrt(Δx^2 + Δy^2)

        # this check isn't really necessary unless user starts the particles in the same place
        if Δr > 0.0
            coef = -hb.k_bond * (1 - hb.l_rest / Δr)
            f_x = coef * Δx
            f_y = coef * Δy

            particle1.f_x += f_x
            particle1.f_y += f_y
            particle2.f_x -= f_x
            particle2.f_y -= f_y
        end
    end
end
