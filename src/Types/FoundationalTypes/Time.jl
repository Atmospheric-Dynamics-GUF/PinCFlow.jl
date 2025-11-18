"""
```julia
Time{A <: Integer, B <: NTuple{3, <:AbstractFloat}}
```

Time integration parameters for the low-storage third-order Runge-Kutta scheme.

```julia
Time()::Time
```

Construct a `Time` instance.

# Fields

  - `nstages::A`: Number of Runge-Kutta stages.

  - `alphark::B`: Runge-Kutta coefficients for the total tendency, i.e. ``\\boldsymbol{\\alpha}_\\mathrm{RK} = \\left(0, - 5 / 9, - 153 / 128\\right)``.

  - `betark::B`: Runge-Kutta coefficients for the previous tendency, i.e. ``\\boldsymbol{\\beta}_\\mathrm{RK} = \\left(1 / 3, 15 / 16, 8 / 15\\right)``.

  - `stepfrac::B`: Time step fractions for each stage, i.e. ``\\boldsymbol{f}_\\mathrm{RK} = \\left(1 / 3, 5 / 12, 1 / 4\\right)``.
"""
struct Time{A <: Integer, B <: NTuple{3, <:AbstractFloat}}
    nstages::A
    alphark::B
    betark::B
    stepfrac::B
end

function Time()::Time
    nstages = 3
    alphark = (0.0, -5.0 / 9.0, -153.0 / 128.0)
    betark = (1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0)
    stepfrac = (1.0 / 3.0, 5.0 / 12.0, 1.0 / 4.0)

    return Time(nstages, alphark, betark, stepfrac)
end
