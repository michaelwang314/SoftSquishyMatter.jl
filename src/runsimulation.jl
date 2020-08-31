export run_simulation!
export save_simulation
export load_simulation

"""
    hr_min_sec(time)

Converts `time` (seconds) to 00:00:00 (hr:min:sec) format.
"""
@inline function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0" : "", hours, 
                  minutes < 10 ? ":0" : ":", minutes, 
                  seconds < 10 ? ":0" : ":", seconds)
end

"""
    run_simulation(simulation; message_interval, save_as)

Run the simulation.  `message_interval` (seconds) controls how often a time
update is printed.  `save_as` is the file to which `simulation` is saved.
"""
function run_simulation!(simulation::Simulation; message_interval::Float64 = 10.0, save_as::String = "", overwrite::Bool = true)
    print_message("SIMULATION STARTED")
    println("Number of threads available: ", Threads.nthreads())
    println("Number of particles: ", length(simulation.particles))
    if length(simulation.bonds) > 0
        println("Number of bonds: ", length(simulation.bonds))
    end
    if length(simulation.angles) > 0
        println("Number of angles: ", length(simulation.angles))
    end
    println("Description: ", simulation.descriptor)
    println("")

    period_x = simulation.periodic_in_x ? simulation.L_x : -1.0
    period_y = simulation.periodic_in_y ? simulation.L_y : -1.0

    if overwrite
        simulation.history = NamedTuple{(:particles, :bonds), Tuple{Array{Particle, 1}, Array{NTuple{2, Particle}, 1}}}[]
        push!(simulation.history, deepcopy(simulation.things_to_save)) # save initial data
    end

    print_message("SIMULATION PROGRESS")
    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    for step = 1 : simulation.num_steps
        for interaction in simulation.interactions
            compute_interactions!(interaction; period_x = period_x, period_y = period_y)
        end
        for external_force in simulation.external_forces
            compute_external_forces!(external_force)
        end

        if step % simulation.save_interval == 0
            push!(simulation.history, deepcopy(simulation.things_to_save))
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
            rate = (step - prev_step) / interval_time
            println(hr_min_sec(time_elapsed), " | ",
                    step, "/", simulation.num_steps, " (", round(step / simulation.num_steps * 100, digits = 1), "%) | ",
                    round(rate, digits = 1), " steps/s | ",
                    hr_min_sec((simulation.num_steps - step) / rate))
            prev_step = step
            interval_start = time()
        end
    end
    println("Average steps/s: ", round(simulation.num_steps / time_elapsed, digits = 1))

    if !isempty(save_as)
        save_simulation(simulation; save_as = save_as)
    end
    print_message("SIMULATION COMPLETED")
end

"""
    save_simulation(simulation; save_as)

Saves `simulation` to `save_as`.
"""
function save_simulation(simulation::Simulation; save_as::String)
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    open(save_as, "w") do f
        serialize(f, simulation)
    end
    println("\nSimulation saved to $save_as\n")
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
    println("\n$file loaded\n")
    return simulation
end
