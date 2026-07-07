"""
```julia
apply_bicgstab!(state::State, tolref::AbstractFloat)::Tuple{Bool, <:Integer}
```

Solve the Poisson equation using a preconditioned BiCGSTAB algorithm and return a tuple containing an error flag and the number of iterations.

# Arguments

  - `state`: Model state.

  - `tolref`: Reference tolerance for convergence criterion.

# See also

  - [`PinCFlow.PoissonSolver.apply_operator!`](@ref)

  - [`PinCFlow.PoissonSolver.apply_preconditioner!`](@ref)

  - [`PinCFlow.MPIOperations.compute_global_dot_product`](@ref)
"""
function apply_bicgstab! end

function apply_bicgstab!(
    state::State,
    tolref::AbstractFloat,
)::Tuple{Bool, <:Integer}
    (; x_size, y_size, z_size) = state.namelists.domain
    (; tolerance, poisson_iterations, preconditioner, tolerance_is_relative) =
        state.namelists.poisson
    (; master, comm, nx, ny, nz) = state.domain
    (; rhs, solution) = state.poisson
    (; p, r0, rold, r, s, t, v, matvec, v_pc) = state.poisson.bicgstab

    # Print information.
    if master
        println(repeat("-", 80))
        println("BiCGSTAB: Solving linear system...")
        println(repeat("-", 80))
        println("")
    end

    # Initialize solution.
    @share solution = 0.0

    # Set parameters.
    maxit = poisson_iterations

    if tolerance_is_relative
        tol = tolerance
    else
        tol = tolerance / tolref
    end

    # Set error flag.
    errflag = false

    apply_operator!(solution, matvec, Total(), state)
    @share for k in 1:nz, j in 1:ny, i in 1:nx
        r0[i, j, k] = rhs[i, j, k] - matvec[i, j, k]
        p[i, j, k] = r0[i, j, k]
        r[i, j, k] = r0[i, j, k]
    end

    res = 0.0
    @share (+, res) for k in 1:nz, j in 1:ny, i in 1:nx
        res += r[i, j, k]^2
    end
    res = MPI.Allreduce(res, +, comm)
    res = sqrt(res / x_size / y_size / z_size)

    b_norm = res

    if res == 0.0 || res / b_norm <= tol
        if master
            println("=> No iteration needed!")
            println("")
        end
        niter = 0
        return (errflag, niter)
    end

    # Loop

    j_b = 0
    while j_b < maxit
        j_b += 1

        # v = A*p
        if preconditioner
            apply_preconditioner!(p, v_pc, state)
        else
            @share v_pc = p
        end
        apply_operator!(v_pc, matvec, Total(), state)
        @share v = matvec

        alpha =
            compute_global_dot_product(r, r0, state) /
            compute_global_dot_product(v, r0, state)
        @share s = r - alpha * v

        # t = A*s
        if preconditioner
            apply_preconditioner!(s, v_pc, state)
        else
            @share v_pc = s
        end
        apply_operator!(v_pc, matvec, Total(), state)
        @share t = matvec

        omega =
            compute_global_dot_product(t, s, state) /
            compute_global_dot_product(t, t, state)
        @share for k in 1:nz, j in 1:ny, i in 1:nx
            solution[i, j, k] += alpha * p[i, j, k] + omega * s[i, j, k]
            rold[i, j, k] = r[i, j, k]
            r[i, j, k] = s[i, j, k] - omega * t[i, j, k]
        end

        #-----------------------
        #   Abort criterion
        #-----------------------

        res = 0.0
        @share (+, res) for k in 1:nz, j in 1:ny, i in 1:nx
            res += r[i, j, k]^2
        end
        res = MPI.Allreduce(res, +, comm)
        res = sqrt(res / x_size / y_size / z_size)

        if res / b_norm <= tol
            if master
                println("Iterations: ", j_b)
                println("Final residual: ", res / b_norm)
                println("")
            end

            niter = j_b

            if preconditioner
                @share s = solution
                apply_preconditioner!(s, solution, state)
            end

            return (errflag, niter)
        end

        beta =
            alpha / omega * compute_global_dot_product(r, r0, state) /
            compute_global_dot_product(rold, r0, state)
        @share p = r + beta * (p - omega * v)
    end

    errflag = true
    niter = j_b

    return (errflag, niter)
end
