export run_simulation!
export save_simulation
export load_simulation

"""
    hr_min_sec(time)

Converts `time` (seconds) to 00:00:00 (hr:min:sec) format.
"""
function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0$hours" : "$hours", ":", minutes < 10 ? "0$minutes" : "$minutes", ":", seconds < 10 ? "0$seconds" : "$seconds")
end

"""
    run_simulation(simulation; message_interval, save_as)

Run the simulation.  `message_interval` (seconds) controls how often a time
update is printed.  `save_as` is the file to which `simulation` is saved.
"""
function run_simulation!(simulation::Simulation; message_interval::Float64 = 10.0, save_as::String = "")
    println("")
    println("   +++++ SIMULATION STARTED +++++")
    println("")

    println("Number of threads available: ", Threads.nthreads())
    println("")
    println("Number of particles: ", length(simulation.particles))
    println("Description: ", simulation.descriptor)
    println("")

    period_x = simulation.periodic_in_x ? simulation.L_x : -1.0
    period_y = simulation.periodic_in_y ? simulation.L_y : -1.0

    simulation.history = Array{Array{Particle, 1}, 1}()

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

        for interaction in simulation.interactions
            compute_interactions!(interaction; period_x = period_x, period_y = period_y)
        end
        for external_force in simulation.external_forces
            compute_external_forces!(external_force)
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
            prev_step = step
            interval_start = time()
        end
    end

    if !isempty(save_as)
        save_simulation(simulation; save_as = save_as)
    end

    println("")
    println("   +++++ SIMULATION COMPLETE +++++")
    println("")
end

"""
    save_simulation(simulation; file)

Saves `simulation` to `file`.
"""
function save_simulation(simulation::Simulation; save_as::String)
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    open(save_as, "w") do f
        serialize(f, simulation)
    end

    println("")
    println("Simulation saved to ", save_as)
    println("")
end

"""
    load_simulation(; file)

Loads simulation from `file`
"""
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
