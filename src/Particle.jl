export group_by_type
export square_lattice_with_φ
export remove_overlaps!

#=
################################################################################
Functions related to 'Particle'
################################################################################
=#

#=
Group particles by Particle.ptype
=#
function group_by_type(particles::Array{Particle, 1}; ptype::Union{Symbol, Array{Symbol, 1}})
    if !(ptype isa Array)
        ptype = [ptype]
    end

    pgroup = Array{Particle, 1}()
    for particle in particles
        if particle.ptype in ptype
            push!(pgroup, particle)
        end
    end
    return pgroup
end

#=
Generate a square lattice of particles with approximate packing fraction
=#
function square_lattice_with_φ(; φ::Float64, L_x::Float64, L_y::Float64, radius::Float64)
    N = L_x * L_y / (pi * radius^2) * φ
    N_x = round(Int64, sqrt(N * L_x / L_y))
    N_y = round(Int64, N_x * L_y / L_x)
    a_x = L_x / N_x
    a_y = L_y / N_y

    positions = Array{Array{Float64, 1}, 1}()
    for i = 0 : N_x - 1
        for j = 0 : N_y - 1
            push!(positions, [i * a_x, j * a_y])
        end
    end
    return positions
end

#=
Remove particles that are within a minimum distance of certain particles
=#
function remove_overlaps!(positions::Array{Array{Float64, 1}, 1}; fixed_positions::Array{Array{Float64, 1}, 1}, minimum_distance::Float64)
    new_positions = Array{Array{Float64, 1}, 1}()
    for (x, y) in positions
        overlap = false
        for (xx, yy) in fixed_positions
            if (x - xx)^2 + (y - yy)^2 < minimum_distance^2
                overlap = true
                break
            end
        end

        if !overlap
            push!(new_positions, [x, y])
        end
    end
    positions = new_positions
end
