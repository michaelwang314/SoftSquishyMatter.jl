using Plots

export plot_frames!
export animate_frames!

"""
    plot_frames!(simulation; frame_nums, frame_size, xlim, ylim, colors, save_to)

Plots each frame of `simulation.history`.
"""
function plot_frames!(simulation::Simulation; frame_nums::Union{Array{Int64, 1}, Nothing} = nothing,
                                              frame_size::Tuple{Int64, Int64} = (600, 600),
                                              xlim::Union{Array{Float64, 1}, Nothing} = nothing, ylim::Union{Array{Float64, 1}, Nothing} = nothing,
                                              colors::Dict{Symbol, String} = Dict(:particle => "black"),
                                              save_to::String = "frames/")
    if !isdir(dirname(save_to))
        mkpath(dirname(save_to))
    end

    frame_nums = isnothing(frame_nums) ? [f for f = 1 : length(simulation.history)] : frame_nums
    xlim = isnothing(xlim) ? [0.0, simulation.L_x] : xlim
    ylim = isnothing(ylim) ? [0.0, simulation.L_y] : ylim

    θ = LinRange(0, 2 * pi, 20)
    unit_circle = (xs = cos.(θ), ys = sin.(θ))

    println("")
    println("   +++++ GENERATING IMAGES +++++")
    println("")

    for f in frame_nums
        plot(size = frame_size, legend = false, axis = false, grid = false, aspect_ratio = 1, xlim = xlim, ylim = ylim)
        for particle in simulation.history[f]
            x, y = particle.x, particle.y
            cx = particle.R * unit_circle.xs .+ x
            cy = particle.R * unit_circle.ys .+ y

            plot!(cx, cy, seriestype = [:shape,], color = colors[particle.ptype], fillalpha = 0.3)
            if !isnothing(particle.active_force)
                af_x, af_y = get_active_force(particle.active_force)
                scale = particle.R / sqrt(af_x^2 + af_y^2)
                plot!([x, x + scale * af_x], [y, y + scale * af_y], color = "red")
            end
        end
        save_as = string(save_to, "Frame $f.png")
        savefig(save_as)
        println(save_as, " done")
    end

    println("")
    println("   +++++ IMAGES GENERATED +++++")
    println("")
end

"""
    animate_frames!(simulation; frame_num, frame_size, xlim, ylim, colors, fps, save_as)

Generates a gif of the simulation.
"""
function animate_frames!(simulation::Simulation; frame_nums::Union{Array{Int64, 1}, Nothing} = nothing,
                                                 frame_size::Tuple{Int64, Int64} = (600, 600),
                                                 xlim::Union{Array{Float64, 1}, Nothing} = nothing, ylim::Union{Array{Float64, 1}, Nothing} = nothing,
                                                 colors::Dict{Symbol, String} = Dict(:particle => "black"),
                                                 fps::Int64 = 10,
                                                 save_as::String = "frames/simulation.gif")
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    frame_nums = isnothing(frame_nums) ? [f for f = 1 : length(simulation.history)] : frame_nums
    xlim = isnothing(xlim) ? [0.0, simulation.L_x] : xlim
    ylim = isnothing(ylim) ? [0.0, simulation.L_y] : ylim

    println("")
    println("   +++++ GENERATING ANIMATION +++++")
    println("")

    θ = LinRange(0, 2 * pi, 20)
    circle = (xs = cos.(θ), ys = sin.(θ))
    animation = @animate for f in frame_nums
        new_frame = true
        for particle in simulation.history[f]
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
    scale = particle.R / sqrt(af_x^2 + af_y^2)
    lx = [particle.x, particle.x + scale * af_x]
    ly = [particle.y, particle.y + scale * af_y]

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
