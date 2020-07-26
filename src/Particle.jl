export group_by_type
export rectangular_lattice
export triangular_lattice
export wedge
export flat_edge
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
Generate rectangular lattice
=#
function rectangular_lattice(; s_x::Float64, s_y::Float64, M_x::Int64, M_y::Int64)
    positions = Array{Array{Float64, 1}, 1}()
    for i = 0 : M_x - 1, j = 0 : M_y - 1
        push!(positions, [i * s_x, j * s_y])
    end
    return positions, M_x * s_x, M_y * s_y
end

#=
Generate triangular lattice
=#
function triangular_lattice(; s::Float64, M_x::Int64, M_y::Int64)
    cell = [[0.0, 0.0], [s / 2, s * sqrt(3) / 2]]
    positions = Array{Array{Float64, 1}, 1}()
    for i = 0 : M_x - 1, j = 0 : M_y - 1
        push!(positions, cell[1] + [i * s, j * s * sqrt(3)])
        push!(positions, cell[2] + [i * s, j * s * sqrt(3)])
    end
    return positions, M_x * s, M_y * s * sqrt(3)
end

#=
Generage a line of particles
=#
function flat_edge(; from::Array{Float64, 1}, to::Array{Float64, 1}, spacing::Float64)
    edge_length = sqrt((to[1] - from[1])^2 + (to[2] - from[2])^2)
    M = ceil(Int64, edge_length / spacing)

    edge = [from]
    for i = 1 : M
        push!(edge, from + i / M * (to - from))
    end
    return edge
end

#=
Generate a funnel with certain angle and length
=#
function wedge(; angle::Float64, orientation::Float64 = 0.0, edge_length::Float64, spacing::Float64, tip::Array{Float64, 1})
    M = ceil(Int64, edge_length / spacing)
    approx_spacing = edge_length / M

    funnel = [tip]
    for i = 1 : M
        push!(funnel, tip + i * approx_spacing * [cos(angle / 2 + orientation), sin(angle / 2 + orientation)])
        push!(funnel, tip + i * approx_spacing * [cos(-angle / 2 + orientation), sin(-angle / 2 + orientation)])
    end
    return funnel
end

#=
Remove particles that are within a minimum distance of certain particles
=#
function remove_overlaps!(positions::Array{Array{Float64, 1}, 1}; fixed_positions::Array{Array{Float64, 1}, 1}, minimum_distance::Float64, period_x::Float64, period_y::Float64)
    indices_to_remove = Array{Int64, 1}()
    for (n, (x, y)) in enumerate(positions)
        for (xx, yy) in fixed_positions
            Δx = x - xx
            Δy = y - yy

            Δx = abs(Δx) > period_x / 2 ? Δx - sign(Δx) * period_x : Δx
            Δy = abs(Δy) > period_y / 2 ? Δy - sign(Δy) * period_y : Δy

            if Δx^2 + Δy^2 < minimum_distance^2
                push!(indices_to_remove, n)
                break
            end
        end
    end
    deleteat!(positions, indices_to_remove)
end
