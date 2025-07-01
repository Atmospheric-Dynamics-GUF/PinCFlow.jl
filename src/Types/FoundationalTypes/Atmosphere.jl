"""
```julia
Atmosphere{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
```

A type representing the atmospheric background state in terrain-following coordinates.

# Type Parameters

  - `A <: AbstractArray{<:AbstractFloat, 3}`: 3D array type for field variables
  - `B <: AbstractVector{<:AbstractFloat}`: Vector type for parameters

# Fields

  - `pstrattfc::A`: Background pressure [p/pref]
  - `thetastrattfc::A`: Background potential temperature [θ/θref]
  - `rhostrattfc::A`: Background density [ρ/ρref]
  - `bvsstrattfc::A`: Brunt-Väisälä frequency squared [N²/Nref²]
  - `fc::B`: Coriolis parameter [f/fref]

## Parameters

  - `namelists::Namelists`: Configuration parameters
  - `constants::Constants`: Physical constants
  - `domain::Domain`: Computational domain information
  - `grid::Grid`: Grid configuration
  - `model`: Atmospheric model type (e.g., `Boussinesq`, `AbstractModel`)
  - `background`: Background state type (e.g., `Isothermal`, `UniformBoussinesq`)

# Model Types

  - `Boussinesq`: Boussinesq approximation
  - Subtypes of `AbstractModel`: Other atmospheric models

# Background Types

  - `Isothermal`: Isothermal background state
  - `UniformBoussinesq`: Uniform background for Boussinesq model
  - `StratifiedBoussinesq`: Stratified background for Boussinesq model

# Examples

```julia
# Create atmosphere with default model and background
atm = Atmosphere(namelists, constants, domain, grid)

# Create Boussinesq atmosphere with uniform background
atm = Atmosphere(
    namelists,
    constants,
    domain,
    grid,
    Boussinesq(),
    UniformBoussinesq(),
)
```

# Notes

  - All fields are non-dimensionalized using reference values from `Constants`
  - Coordinates are in terrain-following system
  - Grid decomposition for parallel computation is handled through `Domain`
"""
struct Atmosphere{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
    pstrattfc::A
    thetastrattfc::A
    rhostrattfc::A
    bvsstrattfc::A
    fc::B
end

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)
```

Creates an `Atmosphere` instance with model and background types specified in namelists.

# Arguments

  - `namelists`: Configuration settings including:

      + `atmosphere`: Background state parameters
      + `domain`: Domain configuration
      + `setting`: Model settings

  - `constants`: Physical constants and reference values
  - `domain`: Grid domain information including parallel decomposition
  - `grid`: Grid metrics and coordinate information
"""
function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)
    (; model) = namelists.setting
    (; background) = namelists.atmosphere
    return Atmosphere(namelists, constants, domain, grid, model, background)
end

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::UniformBoussinesq,
)
```

Creates a uniform Boussinesq atmosphere with:

  - Constant density (ρ = ρref)
  - Uniform potential temperature (θ = θ0)
  - Zero buoyancy frequency (N² = 0)
  - Coriolis parameter based on specified mode
"""
function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::UniformBoussinesq,
)
    (; theta0_dim, coriolis_mode) = namelists.atmosphere
    (; thetaref) = constants
    (; nxx, nyy, nzz) = domain

    # Set the background fields.
    rhostrattfc = ones(nxx, nyy, nzz)
    thetastrattfc = theta0_dim ./ thetaref .* ones(nxx, nyy, nzz)
    pstrattfc = rhostrattfc .* thetastrattfc
    bvsstrattfc = zeros(nxx, nyy, nzz)

    # Set the Coriolis parameter.
    fc = compute_coriolis_parameter(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::StratifiedBoussinesq,
)
```

Creates a stratified Boussinesq atmosphere with:

  - Constant density: ρ = ρ_ref
  - Uniform potential temperature: θ = θ_0
  - Constant pressure from equation of state: p = ρΘ
  - Constant buoyancy frequency (N² = N0²)
  - Coriolis parameter based on specified mode
"""
function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::StratifiedBoussinesq,
)
    (; buoyancy_frequency, theta0_dim, coriolis_mode) = namelists.atmosphere
    (; tref, thetaref) = constants
    (; nxx, nyy, nzz) = domain

    # Set the background fields.
    rhostrattfc = ones(nxx, nyy, nzz)
    thetastrattfc = theta0_dim ./ thetaref .* ones(nxx, nyy, nzz)
    pstrattfc = rhostrattfc .* thetastrattfc
    bvsstrattfc = (buoyancy_frequency .* tref) .^ 2 .* ones(nxx, nyy, nzz)

    # Set the Coriolis parameter.
    fc = compute_coriolis_parameter(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::AbstractModel,
    background::Isothermal,
)
```

Creates an isothermal atmosphere with:

  - Exponential pressure profile: p(z) = p0 * exp(-σz/γT0)
  - Temperature profile: θ(z) = T0 * exp(κσz/T0)
  - Density from equation of state: ρ = p/θ
  - N² computed from vertical θ gradient
  - Handles boundary conditions for N² calculation
"""
function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::AbstractModel,
    background::Isothermal,
)
    (; nbz) = namelists.domain
    (; zboundaries) = namelists.setting
    (; temp0_dim, press0_dim, coriolis_mode) = namelists.atmosphere
    (; thetaref, pref, kappa, sig, gamma, g_ndim) = constants
    (; sizezz, nxx, nyy, nzz, ko, k0, k1) = domain
    (; ztfc, jac, dz) = grid

    # Initialize the background fields.
    (pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc) =
        (zeros(nxx, nyy, nzz) for i in 1:4)

    t0 = temp0_dim / thetaref
    p0 = press0_dim / pref

    # Compute the background fields.
    pstrattfc .= p0 .* exp.(-sig .* ztfc ./ gamma ./ t0)
    thetastrattfc .= t0 .* exp.(kappa .* sig ./ t0 .* ztfc)
    rhostrattfc .= pstrattfc ./ thetastrattfc

    # Compute the squared buoyancy frequency.
    bvsstrattfc .= 0.0
    for k in k0:k1
        bvsstrattfc[:, :, k] .=
            g_ndim ./ thetastrattfc[:, :, k] ./ jac[:, :, k] .* 0.5 .*
            (thetastrattfc[:, :, k + 1] .- thetastrattfc[:, :, k - 1]) ./ dz
    end

    # Compute the squared buoyancy frequency at the boundaries.
    set_vertical_boundaries_of_field!(
        bvsstrattfc,
        namelists,
        domain,
        zboundaries,
        +,
    )
    if ko == 0
        for k in 1:nbz
            bvsstrattfc[:, :, k] .=
                g_ndim ./ thetastrattfc[:, :, k0 - 1] ./ jac[:, :, k0 - 1] .*
                (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
        end
    end
    if ko + nzz == sizezz
        for k in 1:nbz
            bvsstrattfc[:, :, k1 + k] .=
                g_ndim ./ thetastrattfc[:, :, k1 + 1] ./ jac[:, :, k1 + 1] .*
                (thetastrattfc[:, :, k1 + 1] .- thetastrattfc[:, :, k1]) ./ dz
        end
    end

    # Set the Coriolis parameter.
    fc = compute_coriolis_parameter(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end
