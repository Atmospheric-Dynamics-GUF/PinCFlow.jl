function compute_local_array!(namelists::Namelists, domain::Domain)

  # Get all necessary  fields.
  (; nprocx, nprocy) = namelists.domain
  (; comm, master, nx, ny, nz, local_array, master_array, global_array) = domain

  # Scatter.
  if master
    m = 0
    for ip in 1:nprocx, jp in 1:nprocy, k in 1:nz, j in 1:ny, i in 1:nx
      m += 1
      jg = (jp - 1) * ny + j
      ig = (ip - 1) * nx + i
      master_array[m] = global_array[ig, jg, k]
    end
  end
  MPI.Barrier(comm)
  @views MPI.Scatter!(master_array, local_array[:], comm)

  # Return.
  return
end
