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
    θ = LinRange(0, 2 * pi, 20)
    circle = (xs = cos.(θ), ys = sin.(θ))
    for particle in particles
        x = particle.x
        y = particle.y

        cx = particle.R * circle.xs .+ x
        cy = particle.R * circle.ys .+ y

        color = colors[particle.ptype]
        plot!(cx, cy, seriestype = [:shape,], color = color, linecolor = color, fillalpha = 0.3)

        if !isnothing(particle.active_force)
            af_x, af_y = get_active_force(particle.active_force)
            mag = sqrt(af_x^2 + af_y^2)
            plot!([x, x + particle.R * af_x / mag], [y, y + particle.R * af_y / mag], color = "red")
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
function animate_frames!(history::Array{Array{Particle, 1}, 1}; frame_nums::Union{Array{Int64, 1}, Nothing} = nothing, frame_size::Tuple{Int64, Int64} = (600, 600), xlim::Array{Float64, 1}, ylim::Array{Float64, 1}, colors::Dict{Symbol, String} = Dict(:particle => "black"), fps::Int64 = 10, save_as::String = "simulation.gif")
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    if isnothing(frame_nums)
        frame_nums = [f for f = 1 : length(history)]
    end

    println("")
    println("   +++++ GENERATING ANIMATION +++++")
    println("")
    θ = LinRange(0, 2 * pi, 20)
    circle = (xs = cos.(θ), ys = sin.(θ))
    animation = @animate for f in frame_nums
        new_frame = true
        for particle in history[f]
            draw = new_frame ? drawparticle : drawparticle!
            draw(particle, circle, colors[particle.ptype], frame_size, xlim, ylim)
            new_frame = false

            if !isnothing(particle.active_force)
                drawactiveforce!(particle, frame_size, xlim, ylim)
            end
        end

        println("Frame $f done")
    end
    gif(animation, save_as, fps = fps)
    println("")
    println("   +++++ ANIMATION GENERATED +++++")
    println("")
end

@userplot DrawParticle
@recipe function f(info::DrawParticle)
    particle, circle, color, size, xlim, ylim = info.args
    cx = particle.R * circle.xs .+ particle.x
    cy = particle.R * circle.ys .+ particle.y

    linecolor --> color
    linealpha --> 1
    seriestype := :shape
    seriescolor --> color
    seriesalpha --> 0.3

    size --> size
    aspect_ratio --> 1
    label --> false
    axis --> false
    grid --> false
    xlims --> xlim
    ylims --> ylim

    return cx, cy
end

@userplot DrawActiveForce
@recipe function f(info::DrawActiveForce)
    particle, size, xlim, ylim = info.args
    af_x, af_y = get_active_force(particle.active_force)
    mag = sqrt(af_x^2 + af_y^2)
    lx = [particle.x, particle.x + particle.R * af_x / mag]
    ly = [particle.y, particle.y + particle.R * af_y / mag]

    linecolor --> "red"
    linealpha --> 1

    size --> size
    aspect_ratio --> 1
    label --> false
    axis --> false
    grid --> false
    xlims --> xlim
    ylims --> ylim

    return lx, ly
end
