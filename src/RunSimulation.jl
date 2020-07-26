using Serialization

export use_threads
export run_simulation
export save_simulation
export load_simulation

#=
Macro for using multithreading
=#
macro use_threads(multithreaded::Union{Expr, Symbol}, expr::Expr)
    esc(quote
        if $multithreaded
            Threads.@threads $expr
        else
            $expr
        end
    end)
end

#=
Convert seconds to hr:min:sec format
=#
function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0$hours" : "$hours", ":", minutes < 10 ? "0$minutes" : "$minutes", ":", seconds < 10 ? "0$seconds" : "$seconds")
end

#=
Run simulation
=#
function run_simulation(simulation::Simulation; message_interval::Float64 = 10.0, save_to::String = "")
    println("")
    println("   +++++ SIMULATION STARTED +++++")
    println("")
    println("Number of threads available: ", Threads.nthreads())
    println("")
    println("Number of particles: ", length(simulation.particles))
    println("Description: ", simulation.descriptor)
    println("")

    period_x = -1.0
    if simulation.periodic_in_x
        period_x = simulation.L_x
    end
    period_y = -1.0
    if simulation.periodic_in_y
        period_y = simulation.L_y
    end

    if simulation.overwrite
        simulation.history = Array{Array{Particle, 1}, 1}()
    end

    println("")
    println("   +++++ PROGRESS +++++")
    println("")

    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    @time for step = 0 : simulation.num_steps
        if step % simulation.save_interval == 0
            push!(simulation.history, deepcopy(simulation.save_particles))
        end

        for pair_interaction in simulation.pair_interactions
            compute_pair_interaction!(pair_interaction; period_x = period_x, period_y = period_y)
        end
        for external_force in simulation.external_forces
            compute_external_force(external_force)
        end
        for integrator in simulation.integrators
            update_particles!(integrator; period_x = period_x, period_y = period_y)
        end
        for cell_list in simulation.cell_lists
            update_cell_list!(cell_list)
        end

        interval_time = time() - interval_start
        if interval_time > message_interval || step == simulation.num_steps
            time_elapsed += interval_time
            rate = (step + 1 - prev_step) / interval_time
            println(hr_min_sec(time_elapsed), " | ",
                    step + 1, "/", simulation.num_steps + 1, " (", round((step + 1) / simulation.num_steps * 100, digits = 1), "%)", " | ",
                    round(rate, digits = 1), " steps/s", " | ",
                    hr_min_sec((simulation.num_steps - step) / rate))
            interval_start = time()
            prev_step = step
        end
    end

    if !isempty(save_to)
        save_simulation(simulation; file = save_to)
    end

    println("")
    println("   +++++ SIMULATION COMPLETE +++++")
    println("")
end

#=
Save simulation
=#
function save_simulation(simulation::Simulation; file::String)
    if !isdir(dirname(file))
        mkpath(dirname(file))
    end

    open(file, "w") do f
        serialize(f, simulation)
    end
    println("")
    println("Simulation saved to ", file)
    println("")
end

#=
Load Simulation
=#
function load_simulation(; file::String)
    simulation = begin
        open(file, "r") do f
            deserialize(f)
        end
    end
    println("")
    println(file, " loaded")
    println("")
    return simulation
end
