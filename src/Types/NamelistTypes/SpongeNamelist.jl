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
    spongelayer::Bool = false,
    sponge_uv::Bool = false,
    spongeheight::AbstractFloat = 5.0E-1,
    alpharmax::AbstractFloat = 0.0E+0,
    betarmax::AbstractFloat = 1.0E+0,
    lateralsponge::Bool = false,
    spongetype::AbstractSponge = PolynomialSponge(),
    spongeorder::Integer = 1,
    cosmosteps::Integer = 1,
    relax_to_mean::Bool = true,
    perturbation_period::AbstractFloat = 0.0E+0,
    perturbation_amplitude::AbstractFloat = 0.0E+0,
    relaxation_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
)::SpongeNamelist
```

Construct a `SpongeNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `spongelayer::A`: Switch for enabling Rayleigh-damping in the sponges.

  - `sponge_uv::A`: Switch for applying the RHS sponge to the horizontal wind.

  - `spongeheight::B`: Fractional vertical extent of the sponge.

  - `alpharmax::B`: Rayleigh-damping coefficient of the LHS sponge.

  - `betarmax::B`: Rayleigh-damping coefficient of the RHS sponge, multiplied by the time step.

  - `lateralsponge::A`: Switch for the lateral LHS sponge.

  - `spongetype::C`: Profile of the LHS sponge.

  - `spongeorder::D`: Order of the polynomial LHS sponge.

  - `cosmosteps::D`: Factor by which the time step is multiplied in the damping coefficient of the COSMO-like LHS sponge.

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
    spongelayer::A
    sponge_uv::A
    spongeheight::B
    alpharmax::B
    betarmax::B
    lateralsponge::A
    spongetype::C
    spongeorder::D
    cosmosteps::D
    relax_to_mean::A
    perturbation_period::B
    perturbation_amplitude::B
    relaxation_wind::E
end

function SpongeNamelist(;
    spongelayer::Bool = false,
    sponge_uv::Bool = false,
    spongeheight::AbstractFloat = 5.0E-1,
    alpharmax::AbstractFloat = 0.0E+0,
    betarmax::AbstractFloat = 1.0E+0,
    lateralsponge::Bool = false,
    spongetype::AbstractSponge = PolynomialSponge(),
    spongeorder::Integer = 1,
    cosmosteps::Integer = 1,
    relax_to_mean::Bool = true,
    perturbation_period::AbstractFloat = 0.0E+0,
    perturbation_amplitude::AbstractFloat = 0.0E+0,
    relaxation_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
)::SpongeNamelist
    return SpongeNamelist(
        spongelayer,
        sponge_uv,
        spongeheight,
        alpharmax,
        betarmax,
        lateralsponge,
        spongetype,
        spongeorder,
        cosmosteps,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    )
end
