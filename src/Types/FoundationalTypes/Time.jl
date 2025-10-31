"""
```julia
Time{A <: Integer, B <: NTuple{3, <:AbstractFloat}}
```

Time integration parameters for the low-storage third-order Runge-Kutta scheme.

```julia
Time(namelists::Namelists)::Time
```

Construct a `Time` instance.

# Fields

  - `nstages::A`: Number of Runge-Kutta stages.

  - `alphark::B`: Runge-Kutta coefficients for the total tendency, i.e. ``\\boldsymbol{\\alpha}_\\mathrm{RK} = \\left(0, - 5 / 9, - 153 / 128\\right)``.

  - `betark::B`: Runge-Kutta coefficients for the previous tendency, i.e. ``\\boldsymbol{\\beta}_\\mathrm{RK} = \\left(1 / 3, 15 / 16, 8 / 15\\right)``.

  - `stepfrac::B`: Time step fractions for each stage, i.e. ``\\boldsymbol{f}_\\mathrm{RK} = \\left(1 / 3, 5 / 12, 1 / 4\\right)``.

# Arguments

  - `namelists`: Namelists with all model parameters.
"""
struct Time{A <: Integer, B <: NTuple{3, <:AbstractFloat}}
    nstages::A
    alphark::B
    betark::B
    stepfrac::B
end

function Time(namelists::Namelists)::Time
    (; float_type, integer_type) = namelists.discretization

    nstages = integer_type(3)
    alphark = (float_type(0), float_type(-5 / 9), float_type(-153 / 128))
    betark = (float_type(1 / 3), float_type(15 / 16), float_type(8 / 15))
    stepfrac = (float_type(1 / 3), float_type(5 / 12), float_type(1 / 4))

    return Time(nstages, alphark, betark, stepfrac)
end
