"""
```julia
SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractSponge,
    D <: Integer,
    E <: NTuple{3, <:AbstractFloat},
}
```

Namelist for parameters describing the sponge.

```julia
SpongeNamelist(;
    use_sponge::Bool = false,
    damp_horizontal_wind_on_rhs::Bool = false,
    sponge_extent::AbstractFloat = 5.0E-1,
    alpharmax::AbstractFloat = 0.0E+0,
    betarmax::AbstractFloat = 1.0E+0,
    lateral_sponge::Bool = false,
    sponge_type::AbstractSponge = PolynomialSponge(),
    sponge_order::Integer = 1,
    cosmo_steps::Integer = 1,
    relax_to_mean::Bool = true,
    perturbation_period::AbstractFloat = 0.0E+0,
    perturbation_amplitude::AbstractFloat = 0.0E+0,
    relaxation_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
)::SpongeNamelist
```

Construct a `SpongeNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `use_sponge::A`: Switch for enabling Rayleigh-damping in the sponges.

  - `damp_horizontal_wind_on_rhs::A`: Switch for applying the RHS sponge to the horizontal wind.

  - `sponge_extent::B`: Fractional vertical extent of the sponge.

  - `alpharmax::B`: Rayleigh-damping coefficient of the LHS sponge.

  - `betarmax::B`: Rayleigh-damping coefficient of the RHS sponge, multiplied by the time step.

  - `lateral_sponge::A`: Switch for the lateral LHS sponge.

  - `sponge_type::C`: Profile of the LHS sponge.

  - `sponge_order::D`: Order of the polynomial LHS sponge.

  - `cosmo_steps::D`: Factor by which the time step is multiplied in the damping coefficient of the COSMO-like LHS sponge.

  - `relax_to_mean::A`: Switch for relaxing the wind towards its averages on the terrain-following surfaces. If set to `false`, the wind is relaxed towards `relaxation_wind`.

  - `perturbation_period::B`: Period of an oscillating perturbation on top of `relaxation_wind`.

  - `perturbation_amplitude::B`: Amplitude of an oscillating perturbation on top of `relaxation_wind`.

  - `relaxation_wind::E`: Wind to be obtained through Rayleigh damping in the LHS sponge.
"""
struct SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractSponge,
    D <: Integer,
    E <: NTuple{3, <:AbstractFloat},
}
    use_sponge::A
    damp_horizontal_wind_on_rhs::A
    sponge_extent::B
    alpharmax::B
    betarmax::B
    lateral_sponge::A
    sponge_type::C
    sponge_order::D
    cosmo_steps::D
    relax_to_mean::A
    perturbation_period::B
    perturbation_amplitude::B
    relaxation_wind::E
end

function SpongeNamelist(;
    use_sponge::Bool = false,
    damp_horizontal_wind_on_rhs::Bool = false,
    sponge_extent::AbstractFloat = 5.0E-1,
    alpharmax::AbstractFloat = 0.0E+0,
    betarmax::AbstractFloat = 1.0E+0,
    lateral_sponge::Bool = false,
    sponge_type::AbstractSponge = PolynomialSponge(),
    sponge_order::Integer = 1,
    cosmo_steps::Integer = 1,
    relax_to_mean::Bool = true,
    perturbation_period::AbstractFloat = 0.0E+0,
    perturbation_amplitude::AbstractFloat = 0.0E+0,
    relaxation_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
)::SpongeNamelist
    return SpongeNamelist(
        use_sponge,
        damp_horizontal_wind_on_rhs,
        sponge_extent,
        alpharmax,
        betarmax,
        lateral_sponge,
        sponge_type,
        sponge_order,
        cosmo_steps,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    )
end
