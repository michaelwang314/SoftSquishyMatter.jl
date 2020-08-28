export visualize!

"""
    visualize!(simulation; save_as, fps, frame_nums, frame_size, xlim, ylim, particle_colors, multithreaded)

# Arguments
- `simulation::Simulation`: simulation to visualize
- `save_as::Tuple{Symbol, String}` = (:gif, "temp.gif"):
- `fps::Int64 = 10`: frames per second if creating an animation
- `frame_nums::Union{Array{Int64, 1}, Symbol} = :all`: frame numbers to use.
- `frame_size::Tuple{Int64, Int64} = (600, 600)`: size of frame in pixels
- `xlim::Union{Array{Float64, 1}, Symbol} = :all`: frame limits along x-axis
- `ylim::Union{Array{Float64, 1}, Symbol} = :all`: frame limits along y-axis
- `particle_colors::Union{DefaultDict{Symbol, String, String}, Dict{Symbol, String}} = DefaultDict{Symbol, String}("black")`: particle colors
- `multithreaded::Bool = false`: use threads to update frames
"""
function visualize!(simulation::Simulation; save_as::Tuple{Symbol, String} = (:gif, "temp.gif"),
                                            fps::Int64 = 10,
                                            frame_nums::Union{Array{Int64, 1}, Symbol} = :all,
                                            frame_size::Tuple{Int64, Int64} = (600, 600),
                                            xlim::Union{Array{Float64, 1}, Symbol} = :all, ylim::Union{Array{Float64, 1}, Symbol} = :all,
                                            particle_colors::Union{DefaultDict{Symbol, String, String}, Dict{Symbol, String}} = DefaultDict{Symbol, String}("black"),
                                            multithreaded::Bool = false)
    save_type, save_str = save_as
    if !isdir(dirname(save_str))
        mkpath(dirname(save_str))
    end

    if frame_nums == :all
        frame_nums = [f for f = 1 : length(simulation.history)]
    elseif frame_nums == :start
        frame_nums = [1]
    elseif frame_nums == :end
        frame_nums = [length(simulation.history)]
    elseif frame_nums == :startend
        frame_nums = [1, length(simulation.history)]
    end

    if xlim == :all
        xlim = [0.0, simulation.L_x]
    end
    if ylim == :all
        ylim = [0.0, simulation.L_y]
    end

    if save_type == :gif
        print_message("GENERATING GIF")
        animation = Animation()
    elseif save_type == :images
        print_message("GENERATING IMAGES")
    end

    N_particles = length(simulation.things_to_save.particles)
    N_bonds = length(simulation.things_to_save.bonds)
    scene = plot(size = frame_size, legend = false, axis = false, grid = false, aspect_ratio = 1, xlim = xlim, ylim = ylim)
    for (f, frame_num) = enumerate(frame_nums)
        if f == 1
            θ = LinRange(0.0, 2 * pi, 20)
            unit_circle = (xs = cos.(θ), ys = sin.(θ))
            for particle in simulation.history[1].particles
                cxs, cys = particle.R .* unit_circle.xs .+ particle.x, particle.R .* unit_circle.ys .+ particle.y
                color = particle_colors[particle.ptype]
                plot!(cxs, cys, seriestype = [:shape,], color = color, linecolor = color, fillalpha = 0.3)
            end

            for (particle_1, particle_2) in simulation.history[1].bonds
                x_1, y_1 = particle_1.x, particle_1.y
                x_2, y_2 = particle_2.x, particle_2.y
                Δx = wrap_displacement(x_1 - x_2; period = simulation.periodic_in_x ? simulation.L_x : -1.0)
                Δy = wrap_displacement(y_1 - y_2; period = simulation.periodic_in_y ? simulation.L_y : -1.0)
                plot!([x_2, x_2 + Δx], [y_2, y_2 + Δy], color = "black", linealpha = 0.3)
                plot!([x_1, x_1 - Δx], [y_1, y_1 - Δy], color = "black", linealpha = 0.3)
            end
        else
            prev_frame_num = frame_nums[f - 1]
            for n = 1 : N_particles
                prev_particle, particle = simulation.history[prev_frame_num].particles[n], simulation.history[frame_num].particles[n]
                Δx, Δy = particle.x - prev_particle.x, particle.y - prev_particle.y
                scene.series_list[n][:x] .+= Δx
                scene.series_list[n][:y] .+= Δy
            end
            for n = N_particles + 1 : 2 : N_particles + 2 * N_bonds
                particle_1, particle_2 = simulation.history[frame_num].bonds[div(n - N_particles - 1, 2) + 1]
                x_1, y_1 = particle_1.x, particle_1.y
                x_2, y_2 = particle_2.x, particle_2.y
                Δx = wrap_displacement(x_1 - x_2; period = simulation.periodic_in_x ? simulation.L_x : -1.0)
                Δy = wrap_displacement(y_1 - y_2; period = simulation.periodic_in_y ? simulation.L_y : -1.0)
                scene.series_list[n][:x], scene.series_list[n][:y] = [x_2, x_2 + Δx], [y_2, y_2 + Δy]
                scene.series_list[n + 1][:x], scene.series_list[n + 1][:y] = [x_1, x_1 - Δx], [y_1, y_1 - Δy]
            end
        end

        if save_type == :gif
            frame(animation)
        elseif save_type == :images
            savefig(string(save_str, "Frame ", frame_num, ".png"))
        end
        println("Frame ", frame_num, " done")
    end

    if save_type == :gif
        gif(animation, save_str, fps = fps)
        print_message("GIF GENERATED")
    elseif save_type == :images
        print_message("IMAGES GENERATED")
    end
end
