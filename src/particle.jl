export group_by_type
export rectangular_lattice
export triangular_lattice
export wedge
export line
export remove_overlaps!
export remove_shape!
export nearest_neighbor_pairs

"""
    group_by_type(particles; ptype)

Returns a group of particles with type(s) `ptype`.
"""
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

"""
    rectangular_lattice(; s_x, s_y, M_x, M_y)

Generates a rectangular lattice with lattice constants `s_x, s_y` in the x, y
directions.  Each cell is duplicated `M_x, M_y` times in the x, y directions.
"""
function rectangular_lattice(; s_x::Float64, s_y::Float64, M_x::Int64, M_y::Int64)
    positions = Array{Array{Float64, 1}, 1}()
    for i = 0 : M_x - 1, j = 0 : M_y - 1
        push!(positions, [i * s_x, j * s_y])
    end
    return positions, M_x * s_x, M_y * s_y
end

"""
    triangular_lattice(; s, M_x, M_y)

Generates a triangular lattice with lattice constant `s`.  Each cell consists of
points `(0, 0)` and `(s / 2, s * sqrt(3) / 2)`.  The cells are duplicated
`M_x, M_y` times in the x, y directions.  Note that there are other
differently oriented triangular lattices.
"""
function triangular_lattice(; s::Float64, M_x::Int64, M_y::Int64)
    cell = [[0.0, 0.0], [s / 2, s * sqrt(3) / 2]]
    positions = Array{Array{Float64, 1}, 1}()
    for i = 0 : M_x - 1, j = 0 : M_y - 1
        push!(positions, cell[1] + [i * s, j * s * sqrt(3)])
        push!(positions, cell[2] + [i * s, j * s * sqrt(3)])
    end
    return positions, M_x * s, M_y * s * sqrt(3)
end

"""
    line(; from, to, spacing)

Generates a line of points from point `from` to point `to` with approximate
spacing `spacing`.  Note that the spacing is automatically adjusted slightly so
that points are equally spaced and include the start and end points.
"""
function line(; from::Array{Float64, 1}, to::Array{Float64, 1}, spacing::Float64)
    edge_length = sqrt((to[1] - from[1])^2 + (to[2] - from[2])^2)
    M = ceil(Int64, edge_length / spacing)

    edge = [from]
    for i = 1 : M
        push!(edge, from + i / M * (to - from))
    end
    return edge
end

"""
    wedge(; angle, edge_length, spacing, orientation, tip)

Generates a V-shaped wedge.  The angle of the wedge is `angle` with each arm a
length `edge_length`.  The points are separated by approximately `spacing` so
that the start and end points are included.  `orientation` is the angle relative
to the x-axis and tip is the x, y location of the tip of the wedge.
"""
function wedge(; angle::Float64, edge_length::Float64, spacing::Float64, orientation::Float64 = 0.0, tip::Array{Float64, 1})
    M = ceil(Int64, edge_length / spacing)
    approx_spacing = edge_length / M

    funnel = [tip]
    for i = 1 : M
        push!(funnel, tip + i * approx_spacing * [cos(angle / 2 + orientation), sin(angle / 2 + orientation)])
        push!(funnel, tip + i * approx_spacing * [cos(-angle / 2 + orientation), sin(-angle / 2 + orientation)])
    end
    return funnel
end

"""
    remove_overlaps!(positions; fixed_positions, minimum_distance, period_x, period_y)

Removes points from `positions` that are within a `minimum_distance` of points
in `fixed_positions`.  This is used to avoid large forces due to overlaps.
"""
function remove_overlaps!(positions::Array{Array{Float64, 1}, 1}; fixed_positions::Array{Array{Float64, 1}, 1}, minimum_distance::Float64, period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    indices_to_remove = Array{Int64, 1}()
    for (n, (x, y)) in enumerate(positions)
        for (xx, yy) in fixed_positions
            Δx = wrap_displacement(x - xx; period = period_x)
            Δy = wrap_displacement(y - yy; period = period_y)

            if Δx^2 + Δy^2 < minimum_distance^2
                push!(indices_to_remove, n)
                break
            end
        end
    end
    deleteat!(positions, indices_to_remove)
end

"""
    remove_shape!(positions; shape, in_or_out)

Remove positions that are either inside (:in) or outside (:out) of a 
parametric shape.
"""
function remove_shape!(positions::Array{Array{Float64, 1}, 1}; shape::Array{Array{Float64, 1}, 1}, in_or_out::Symbol)
    indices_to_remove = Array{Int64, 1}()
    N_shape = length(shape)
    for (n, (x, y)) in enumerate(positions)
        θ = 0.0
        for p = 1 : N_shape
            next_p = p % N_shape + 1
            Δx_1, Δy_1 = shape[p][1] - x, shape[p][2] - y
            Δx_2, Δy_2 = shape[next_p][1] - x, shape[next_p][2] - y
            cross = (Δx_1 * Δy_2 - Δx_2 * Δy_1) / sqrt((Δx_1^2 + Δy_1^2) * (Δx_2^2 + Δy_2^2))

            θ += asin(cross)
        end

        winding_number = round(θ / (2 * pi))
        if (winding_number == 0.0 && in_or_out == :out) || (winding_number != 0.0 && in_or_out == :in)
            push!(indices_to_remove, n)
        end
    end
    deleteat!(positions, indices_to_remove)
end

"""
    nearest_neighbor_pairs(particles; maximum_distance)

Returns pairs of particles that are within `maximum_distance` of eachother.
"""
function nearest_neighbor_pairs(particles::Array{Particle, 1}; maximum_distance::Float64 = 0.0, period_x::Float64 = -1.0, period_y::Float64 = -1.0)
    pairs = Array{Tuple{Particle, Particle}, 1}()

    N = length(particles)
    for i = 1 : N, j = i + 1 : N
        Δx = wrap_displacement(particles[i].x - particles[j].x; period = period_x)
        Δy = wrap_displacement(particles[i].y - particles[j].y; period = period_y)

        if Δx^2 + Δy^2 < maximum_distance^2
            push!(pairs, (particles[i], particles[j]))
        end
    end

    return pairs
end
