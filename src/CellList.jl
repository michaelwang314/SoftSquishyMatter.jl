export update_cell_list!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
    update_cell_list!(cell_list)

Updates cell list with current positions of particles stored under
cell_list.particles.
"""
function update_cell_list!(cell_list::CellList)
    fill!(cell_list.start_pid, -1)

    for (n, particle) in enumerate(cell_list.particles)
        i = trunc(Int64, particle.x / cell_list.cell_spacing_x) + 1
        j = trunc(Int64, particle.y / cell_list.cell_spacing_y) + 1

        if cell_list.start_pid[i, j] > 0
            cell_list.next_pid[n] = cell_list.start_pid[i, j]
        end
        cell_list.start_pid[i, j] = n
    end
end
