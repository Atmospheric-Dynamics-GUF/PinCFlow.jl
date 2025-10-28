"""
```julia
SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: Bool,
    F <: AbstractSponge,
    G <: Integer,
    H <: Integer,
    I <: Bool,
    J <: AbstractFloat,
    K <: AbstractFloat,
    L <: NTuple{3, <:AbstractFloat},
}
```

Namelist for parameters describing the sponge.

```julia
SpongeNamelist(;
    damp_horizontal_wind_on_rhs::Bool = false,
    sponge_extent::AbstractFloat = 5.0E-1,
    alpharmax::AbstractFloat = 0.0E+0,
    betarmax::AbstractFloat = 0.0E+0,
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

  - `damp_horizontal_wind_on_rhs::A`: Switch for applying the RHS sponge to the horizontal wind.

  - `sponge_extent::B`: Fractional vertical extent of the sponge.

  - `alpharmax::C`: Rayleigh-damping coefficient of the LHS sponge.

  - `betarmax::D`: Rayleigh-damping coefficient of the RHS sponge, multiplied by the time step.

  - `lateral_sponge::E`: Switch for the lateral LHS sponge.

  - `sponge_type::F`: Profile of the LHS sponge.

  - `sponge_order::G`: Order of the polynomial LHS sponge.

  - `cosmo_steps::H`: Factor by which the time step is multiplied in the damping coefficient of the COSMO-like LHS sponge.

  - `relax_to_mean::I`: Switch for relaxing the wind towards its averages on the terrain-following surfaces. If set to `false`, the wind is relaxed towards `relaxation_wind`.

  - `perturbation_period::J`: Period of an oscillating perturbation on top of `relaxation_wind`.

  - `perturbation_amplitude::K`: Amplitude of an oscillating perturbation on top of `relaxation_wind`.

  - `relaxation_wind::L`: Wind to be obtained through Rayleigh damping in the LHS sponge.
"""
struct SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: Bool,
    F <: AbstractSponge,
    G <: Integer,
    H <: Integer,
    I <: Bool,
    J <: AbstractFloat,
    K <: AbstractFloat,
    L <: NTuple{3, <:AbstractFloat},
}
    damp_horizontal_wind_on_rhs::A
    sponge_extent::B
    alpharmax::C
    betarmax::D
    lateral_sponge::E
    sponge_type::F
    sponge_order::G
    cosmo_steps::H
    relax_to_mean::I
    perturbation_period::J
    perturbation_amplitude::K
    relaxation_wind::L
end

function SpongeNamelist(;
    damp_horizontal_wind_on_rhs::Bool = false,
    sponge_extent::AbstractFloat = 5.0E-1,
    alpharmax::AbstractFloat = 0.0E+0,
    betarmax::AbstractFloat = 0.0E+0,
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
