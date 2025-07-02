"""
```julia
Time{A <: Integer, B <: AbstractVector{<:AbstractFloat}}
```

Time integration parameters for explicit Runge-Kutta scheme.

This struct encapsulates the coefficients and stage information for multi-stage explicit Runge-Kutta time integration methods used to advance the governing equations in time.

# Type Parameters

  - `A<:Integer`: Integer type for number of stages
  - `B<:AbstractVector{<:AbstractFloat}`: Vector type for RK coefficients

# Fields

  - `nstages::A`: Number of Runge-Kutta stages
  - `alphark::B`: RK coefficients for intermediate stage combinations [α₁, α₂, α₃]
  - `betark::B`: RK coefficients for time step fractions [β₁, β₂, β₃]
  - `stepfrac::B`: Fractional time steps for each stage
"""
struct Time{A <: Integer, B <: AbstractVector{<:AbstractFloat}}
    nstages::A
    alphark::B
    betark::B
    stepfrac::B
end

"""
```julia
Time()
```

Constructs a `Time` instance.
"""
function Time()

    # Set Runge-Kutta parameters.
    nstages = 3
    alphark = [0.0, -5.0 / 9.0, -153.0 / 128.0]
    betark = [1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0]
    stepfrac = [1.0 / 3.0, 5.0 / 12.0, 1.0 / 4.0]

    # Return a Time instance.
    return Time(nstages, alphark, betark, stepfrac)
end
