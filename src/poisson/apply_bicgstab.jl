function apply_bicgstab!(
  b_in::AbstractArray{<:AbstractFloat, 3},
  tolref::AbstractFloat,
  dt::AbstractFloat,
  sol::AbstractArray{<:AbstractFloat, 3},
  res::AbstractFloat,
  niter::Integer,
  errflag::Bool,
  namelists::Namelists,
  domain::Domain,
  grid::Grid,
  poisson::Poisson,
)
  (; sizex, sizey, sizez) = namelists.domain
  (; tolpoisson, maxiterpoisson, preconditioner, relative_tolerance) =
    namelists.poisson
  (; master, comm, nx, ny, nz) = domain
  (; r_vm, b_vm, p, r0, rold, r, s, t, v, matvec, v_pc) = poisson.bicgstab

  # Print information.
  if master
    println("")
    println(repeat("-", 80))
    println("BICGSTAB: Solving linear system...")
    println(repeat("-", 80))
    println("")
  end

  # Initialize solution.
  sol .= 0.0

  # Set parameters.
  maxit = maxiterpoisson

  if relative_tolerance
    tol = tolpoisson
  else
    tol = tolpoisson / tolref
  end

  # Set error flag.
  errflag = false

  b = b_in

  apply_operator!(sol, matvec, Total(), namelists, domain, poisson)
  r0 .= b .- matvec
  p .= r0
  r .= r0

  res_local = 0.0
  for k in 1:nz
    for j in 1:ny
      for i in 1:nx
        res_local = res_local + r[i, j, k]^2
      end
    end
  end

  # MPI: Find global residual.
  res = MPI.Allreduce(res_local, +, comm)

  res = sqrt(res / sizex / sizey / sizez)

  b_norm = res

  r_vm .= 0.0
  for k in 1:nz
    r_vm .= r_vm .+ r[:, :, k]
  end
  r_vm .= r_vm ./ sizez

  res_local = 0.0
  for j in 1:ny
    for i in 1:nx
      res_local = res_local + r_vm[i, j]^2
    end
  end

  res_vm = MPI.Allreduce(res_local, +, comm)

  res_vm = sqrt(res_vm / sizex / sizey)

  b_vm_norm = res_vm

  if res == 0.0 || res / b_norm <= tol
    if master
      println("==> No iteration needed!")
    end
    niter = 0
    return
  end

  # Loop

  for j_b in 1:maxit

    # v = A*p
    if preconditioner
      apply_preconditioner!(p, v_pc, namelists, domain, grid, poisson)
    else
      v_pc .= p
    end
    apply_operator!(v_pc, matvec, Total(), namelists, domain, poisson)
    v .= matvec

    alpha =
      compute_global_dot_product(r, r0) / compute_global_dot_product(v, r0)
    s .= r .- alpha .* v

    # t = A*s
    if preconditioner
      apply_preconditioner!(s, v_pc, namelists, domain, grid, poisson)
    else
      v_pc .= s
    end
    apply_operator!(v_pc, matvec, Total(), namelists, domain, poisson)
    t .= matvec

    omega = compute_global_dot_product(t, s) / compute_global_dot_product(t, t)
    sol .= sol .+ alpha .* p .+ omega .* s

    rold .= r
    r .= s .- omega .* t

    #-----------------------
    #   Abort criterion
    #-----------------------

    res_local = 0.0
    for k in 1:nz
      for j in 1:ny
        for i in 1:nx
          res_local = res_local + r[i, j, k]^2
        end
      end
    end

    # MPI: Find global residual.
    res = MPI.Allreduce(res_local, +, comm)

    res = sqrt(res / sizex / sizey / sizez)

    r_vm .= 0.0
    for k in 1:nz
      r_vm .= r_vm .+ r[:, :, k]
    end
    r_vm .= r_vm ./ sizez

    res_local = 0.0
    for j in 1:ny
      for i in 1:nx
        res_local = res_local + r_vm[i, j]^2
      end
    end

    res_vm = MPI.Allreduce(res_local, +, comm)

    res_vm = sqrt(res_vm / sizex / sizey)

    if max(res / b_norm, res_vm / b_vm_norm) <= tol
      if master
        println("Nb.of iterations: j = ", j_b)
        println("Final residual: res = ", res / b_norm)
        println("Final residual v.m. = ", res_vm / b_vm_norm)
        println("")
      end

      niter = j_b

      if preconditioner
        s .= sol
        apply_preconditioner!(s, sol, namelists, domain, grid, poisson)
      end

      return
    end

    beta =
      alpha / omega * compute_global_dot_product(r, r0) /
      compute_global_dot_product(rold, r0)
    p .= r .+ beta .* (p .- omega .* v)
  end

  # max iteration

  errflag = true
  niter = j_b

  # Return.
  return
end
