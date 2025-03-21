function compute_global_array!(namelists::Namelists, domain::Domain)

  # Get all necessary  fields.
  (; nprocx, nprocy) = namelists.domain
  (; comm, master, nx, ny, nz, local_array, master_array, global_array) = domain

  # Gather.
  for k in 1:nz
    for j in 1:ny
      MPI.Gather!(view(local_array, :, j, k), view(master_array, :, j, k), comm)
    end
    MPI.Barrier(comm)
    if master
      for j in 1:ny
        jm = j
        for jp in 1:nprocy
          jg = ny * (jp - 1) + j
          for ip in 1:nprocx
            for i in 1:nx
              ig = nx * (ip - 1) + i
              im = nprocy * nx * (ip - 1) + (jp - 1) * nx + i
              global_array[ig, jg, k] = master_array[im, jm, k]
            end
          end
        end
      end
    end
  end

  # Return.
  return
end
