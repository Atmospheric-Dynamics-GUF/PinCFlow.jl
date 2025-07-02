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
    SpongeNamelist(; <keyword arguments>)

Configuration parameters for damping layers that absorb waves near domain boundaries.

# Arguments

  - `spongelayer::Bool = false`: Enable sponge layer damping
  - `sponge_uv::Bool = false`: Apply sponge to horizontal velocity components
  - `spongeheight::AbstractFloat = 5.0E-1`: Fractional height of vertical sponge region
  - `spongealphaz_dim::AbstractFloat = 1.0E-2`: Dimensional damping coefficient [1/s]
  - `spongealphaz_fac::AbstractFloat = 1.0E+0`: Damping coefficient scaling factor
  - `unifiedsponge::Bool = false`: Use unified damping formulation for all variables
  - `lateralsponge::Bool = false`: Enable horizontal boundary damping
  - `spongetype::AbstractSponge = PolynomialSponge()`: Spatial damping profile. Options: `ExponentialSponge()`, `COSMOSponge()`, `PolynomialSponge()`, `SinusoidalSponge()`
  - `spongeorder::Integer = 1`: Polynomial order for `PolynomialSponge()`
  - `cosmosteps::Integer = 1`: Time substeps for COSMO-type damping
  - `relax_to_mean::Bool = true`: Damp toward horizontal mean instead of zero
  - `perturbation_period::AbstractFloat = 0.0E+0`: Period for time-varying relaxation [s]
  - `perturbation_amplitude::AbstractFloat = 0.0E+0`: Amplitude of relaxation perturbation
  - `relaxation_wind::NTuple{3,AbstractFloat} = (0,0,0)`: Target wind field [m/s] when not relaxing to mean

# Usage
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
