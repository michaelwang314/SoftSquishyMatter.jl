using Plots

export plot_history!
export plot_frame!

#=
Plot a frame
=#
function plot_frame!(particles::Array{Particle, 1}; frame_size::Tuple{Int64, Int64},
                                                    xlim::Array{Float64, 1}, ylim::Array{Float64, 1},
                                                    colors::Dict{Symbol, String} = Dict(:particle => "black"),
                                                    save_as::String)
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    plot(size = frame_size, legend = false, axis = false, grid = false, aspectratio = 1, xlim = xlim, ylim = ylim)
    θ = LinRange(0, 2 * pi, 20)
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
                plot!([x, x + 1.5 * particle.R * cos(particle.active_force.θ::Float64)], [y, y + 1.5 * particle.R * sin(particle.active_force.θ::Float64)], color = "red")
            end
        end
    end
    savefig(save_as)
    println(save_as, " done")
end

#=
Plot history
=#
function plot_history!(history::Array{Array{Particle, 1}, 1}; frame_size::Tuple{Int64, Int64}, xlim::Array{Float64, 1}, ylim::Array{Float64, 1}, colors::Dict{Symbol, String} = Dict(:particle => "black"), folder::String = "frames/")
    println("")
    println("   +++++ GENERATING IMAGES +++++")
    println("")
    for (f, particles) in enumerate(history[frame_nums])
        plot_frame!(particles; frame_size = frame_size, xlim = xlim, ylim = ylim, colors = colors, save_as = string(folder, "Frame $f.png"))
    end
    println("")
    println("   +++++ IMAGES GENERATED +++++")
    println("")
end
