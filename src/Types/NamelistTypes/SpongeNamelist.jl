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

Namelist for sponge configuration (see constructor for parameter descriptions).
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
    spongealphaz_dim::B
    spongealphaz_fac::B
    unifiedsponge::A
    lateralsponge::A
    spongetype::C
    spongeorder::D
    cosmosteps::D
    relax_to_mean::A
    perturbation_period::B
    perturbation_amplitude::B
    relaxation_wind::E
end

"""
```julia
SpongeNamelist(;
    spongelayer = false,
    sponge_uv = false,
    spongeheight = 5.0E-1,
    spongealphaz_dim = 1.0E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = false,
    lateralsponge = false,
    spongetype = PolynomialSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = true,
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
    relaxation_wind = (0.0E+0, 0.0E+0, 0.0E+0),
)
```

Configuration parameters for damping layers that absorb waves near domain boundaries.

# Arguments

  - `spongelayer`: Enable sponge layer damping
  - `sponge_uv`: Apply sponge to horizontal velocity components
  - `spongeheight:AbstractFloat = `: Fractional height of vertical sponge region
  - `spongealphaz_dim:AbstractFloat = `: Dimensional damping coefficient [1/s]
  - `spongealphaz_fac:AbstractFloat = `: Damping coefficient scaling factor
  - `unifiedsponge`: Use unified damping formulation for all variables
  - `lateralsponge`: Enable horizontal boundary damping
  - `spongetype:AbstractSponge`: Spatial damping profile. Options: `ExponentialSponge()`, `COSMOSponge()`, `PolynomialSponge()`, `SinusoidalSponge()`
  - `spongeorde`: Polynomial order for `PolynomialSponge()`
  - `cosmostep`: Time substeps for COSMO-type damping
  - `relax_to_mean`: Damp toward horizontal mean instead of zero
  - `perturbation_period`: Period for time-varying relaxation [s]
  - `perturbation_amplitude`: Amplitude of relaxation perturbation
  - `relaxation_wind`: Target wind field [m/s] when not relaxing to mean
"""
function SpongeNamelist(;
    spongelayer = false,
    sponge_uv = false,
    spongeheight = 5.0E-1,
    spongealphaz_dim = 1.0E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = false,
    lateralsponge = false,
    spongetype = PolynomialSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = true,
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
    relaxation_wind = (0.0E+0, 0.0E+0, 0.0E+0),
)
    return SpongeNamelist(
        spongelayer,
        sponge_uv,
        spongeheight,
        spongealphaz_dim,
        spongealphaz_fac,
        unifiedsponge,
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
