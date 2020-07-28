using Plots

export plot_frame!
export plot_frames!
export animate_frames!

"""
    plot_frame!(particles; frame_size, xlim, ylim, colors, save_as)

Displays current configuration of `particles`.
"""
function plot_frame!(particles::Array{Particle, 1}; frame_size::Tuple{Int64, Int64} = (600, 600),
                                                    xlim::Array{Float64, 1}, ylim::Array{Float64, 1},
                                                    colors::Dict{Symbol, String} = Dict(:particle => "black"),
                                                    save_as::String)
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    plot(size = frame_size, legend = false, axis = false, grid = false, aspectratio = 1, xlim = xlim, ylim = ylim)
    θ = LinRange(0, 2 * pi, 15)
    circle = [cos.(θ), sin.(θ)]
    for particle in particles
        x = particle.x
        y = particle.y

        cx = particle.R * circle[1] .+ x
        cy = particle.R * circle[2] .+ y

        color = colors[particle.ptype]
        plot!(cx, cy, seriestype = [:shape,], color = color, linecolor = color, fillalpha = 0.3)

        if !isnothing(particle.active_force)
            if (particle.active_force isa ActiveBrownian) || (particle.active_force isa RunAndTumble)
                plot!([x, x + particle.R * cos(particle.active_force.θ::Float64)], [y, y + particle.R * sin(particle.active_force.θ::Float64)], color = "red")
            end
        end
    end
    savefig(save_as)
    println(save_as, " done")
end

"""
    plot_frames!(history; frame_nums, frame_size, xlim, ylim, colors, folder)

Plots each frame of `history`.
"""
function plot_frames!(history::Array{Array{Particle, 1}, 1}; frame_nums::Union{Array{Int64, 1}, Nothing} = nothing, frame_size::Tuple{Int64, Int64} = (600, 600), xlim::Array{Float64, 1}, ylim::Array{Float64, 1}, colors::Dict{Symbol, String} = Dict(:particle => "black"), folder::String = "frames/")
    if isnothing(frame_nums)
        frame_nums = [f for f = 1 : length(history)]
    end

    println("")
    println("   +++++ GENERATING IMAGES +++++")
    println("")
    for f in frame_nums
        plot_frame!(history[f]; frame_size = frame_size, xlim = xlim, ylim = ylim, colors = colors, save_as = string(folder, "Frame $f.png"))
    end
    println("")
    println("   +++++ IMAGES GENERATED +++++")
    println("")
end

"""
    animate_frames!(...)

...
"""
function animate_frames!(history::Array{Array{Particle, 1}, 1}; frame_nums::Union{Array{Int64, 1}, Nothing} = nothing, frame_size::Tuple{Int64, Int64} = (600, 600), xlim::Array{Float64, 1}, ylim::Array{Float64, 1}, colors::Dict{Symbol, String} = Dict(:particle => "black"), save_as::String = "simulation.gif")
    
end
