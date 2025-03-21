function compute_local_array!(namelists::Namelists, domain::Domain)

  # Get all necessary  fields.
  (; nprocx, nprocy) = namelists.domain
  (; comm, master, nx, ny, nz, local_array, master_array, global_array) = domain

  # Scatter.
  for k in 1:nz
    if master
      for j in 1:ny
        jm = j
        for jp in 1:nprocy
          jg = ny * (jp - 1) + j
          for ip in 1:nprocx, i in 1:nx
            ig = nx * (i_prc - 1) + i
            im = nprocy * nx * (ip - 1) + (jp - 1) * nx + i
            master_array[im, jm, k] = global_array[ig, jg, k]
          end
        end
      end
    end
    MPI.Barrier(comm)
    for j in 1:ny
      MPI.Scatter!(
        view(master_array, :, j, k),
        view(local_array, :, j, k),
        comm,
      )
    end
  end

  # Return.
  return
end
