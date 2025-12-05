"""
```julia
WKB{
    A <: Integer,
    B <: AbstractArray{<:Integer, 3},
    C <: Rays,
    D <: MergedRays,
    E <: SurfaceIndices,
    F <: WKBIncrements,
    G <: WKBIntegrals,
    H <: WKBTendencies,
    I <: Ref{<:AbstractFloat},
    J <: AbstractArray{<:AbstractFloat, 3},
    K <: AbstractMatrix{<:AbstractFloat},
}
```

Main container for WKB ray-tracing data and parameters.

```julia
WKB(namelists::Namelists, domain::Domain)::WKB
```

Construct a `WKB` instance by dispatching to a test-case-specific method.

```julia
WKB(namelists::Namelists, domain::Domain, wkb_mode::NoWKB)::WKB
```

Construct a `WKB` instance with zero-size arrays for non-WKB configurations.

```julia
WKB(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::WKB
```

Construct a `WKB` instance.

This method primarily determines the size of the spectral dimension of ray-volume arrays and initializes them and related arrays (with zeros) accordingly. The proper initialization with nonzero wave action is performed by [`PinCFlow.MSGWaM.RayUpdate.initialize_rays!`](@ref).

# Fields

  - `nxray::A`: Number of ray volumes allowed in ``\\widehat{x}``, per grid cell and wave mode (`multiplication_factor * nrx * nrk`, taken from `namelists.wkb`).

  - `nyray::A`: Number of ray volumes allowed in ``\\widehat{y}``, per grid cell and wave mode (`multiplication_factor * nry * nrl`, taken from `namelists.wkb`).

  - `nzray::A`: Number of ray volumes allowed in ``\\widehat{z}``, per grid cell and wave mode (`multiplication_factor * nrz * nrm`, taken from `namelists.wkb`).

  - `nxray_wrk::A`: `2 * nxray`.

  - `nyray_wrk::A`: `2 * nyray`.

  - `nzray_wrk::A`: `2 * nzray`.

  - `nray_max::A`: Maximum ray-volume count allowed per grid-cell before merging is triggered (`nxray * nyray * nzray * namelists.wkb.wave_modes`).

  - `nray_wrk::A`: Size of the spectral dimension of ray-volume arrays (`nxray_wrk * nyray_wrk * nzray_wrk`).

  - `n_sfc::A`: Number of orographic wave modes.

  - `nray::B`: Ray-volume count in each grid cell.

  - `rays::C`: Prognostic ray-volume properties.

  - `merged_rays::D`: Container used for creating merged ray volumes.

  - `surface_indices::E`: Indices that connect orographic wave modes to ray volumes.

  - `increments::F`: WKBIncrements of the prognostic ray-volume properties.

  - `integrals::G`: Integrals of ray-volume properties.

  - `tendencies::H`: Gravity-wave drag and heating fields.

  - `cgx_max::I`: Maximum zonal group velocities.

  - `cgy_max::I`: Maximum meridional group velocities.

  - `cgz_max::J`: Maximum vertical group velocities.

  - `zb::K`: Upper edge of the blocked layer.

  - `diffusion::J`: Diffusion induced by wave breaking.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Types.WKBTypes.Rays`](@ref)

  - [`PinCFlow.Types.WKBTypes.SurfaceIndices`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBIncrements`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBIntegrals`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBTendencies`](@ref)
"""
struct WKB{
    A <: Integer,
    B <: AbstractArray{<:Integer, 3},
    C <: Rays,
    D <: MergedRays,
    E <: SurfaceIndices,
    F <: WKBIncrements,
    G <: WKBIntegrals,
    H <: WKBTendencies,
    I <: Ref{<:AbstractFloat},
    J <: AbstractArray{<:AbstractFloat, 3},
    K <: AbstractMatrix{<:AbstractFloat},
    L <: TriadTendencies
}
    nxray::A
    nyray::A
    nzray::A
    nxray_wrk::A
    nyray_wrk::A
    nzray_wrk::A
    nray_max::A
    nray_wrk::A
    n_sfc::A
    nray::B
    rays::C
    merged_rays::D
    surface_indices::E
    increments::F
    integrals::G
    tendencies::H
    cgx_max::I
    cgy_max::I
    cgz_max::J
    zb::K
    diffusion::J
    spec_tend::L
end

function WKB(namelists::Namelists, domain::Domain)::WKB
    (; wkb_mode) = namelists.wkb

    return WKB(namelists, domain, wkb_mode)
end

function WKB(namelists::Namelists, domain::Domain, wkb_mode::NoWKB)::WKB
    return WKB(
        [0 for i in 1:9]...,
        zeros(Int, 0, 0, 0),
        Rays(0, 0, 0, 0),
        MergedRays(0, 0),
        SurfaceIndices(0, 0, 0),
        WKBIncrements(0, 0, 0, 0),
        WKBIntegrals(0, 0, 0),
        WKBTendencies(0, 0, 0),
        [Ref(0.0) for i in 1:2]...,
        zeros(0, 0, 0),
        zeros(0, 0),
        zeros(0, 0, 0),
        TriadTendencies(0, 0, 0, 0, 0),
    )
end

function WKB(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::WKB
    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        multiplication_factor,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
    ) = namelists.wkb
    (; x_size, y_size, z_size) = namelists.domain
    (; kp_size, m_size) = namelists.triad
    (; nxx, nyy, nzz) = domain

    # Check if spectral-extent factors are set correctly.
    if x_size > 1 && dkr_factor == 0.0
        error("Error in WKB: x_size > 1 && dkr_factor == 0!")
    end
    if y_size > 1 && dlr_factor == 0.0
        error("Error in WKB: y_size > 1 && dlr_factor == 0!")
    end
    if z_size == 1 || dmr_factor == 0.0
        error("Error in WKB: z_size == 1 || dmr_factor == 0!")
    end

    # Set zonal ray-volume count.
    if x_size == 1
        nxray = 1
    else
        nxray = multiplication_factor * nrx * nrk
    end

    # Set meridional ray-volume count.
    if y_size == 1
        nyray = 1
    else
        nyray = multiplication_factor * nry * nrl
    end

    # Set vertical ray-volume count.
    nzray = multiplication_factor * nrz * nrm

    # Set maximum ray-volume count.
    nray_max = nxray * nyray * nzray * wave_modes

    # Set spectral dimension of ray-volume array.
    if nxray > 1
        nxray_wrk = 2 * nxray
    else
        nxray_wrk = 1
    end
    if nyray > 1
        nyray_wrk = 2 * nyray
    else
        nyray_wrk = 1
    end
    if nzray > 1
        nzray_wrk = 2 * nzray
    else
        nzray_wrk = 1
    end
    nray_wrk = nxray_wrk * nyray_wrk * nzray_wrk * wave_modes

    # Set number of surface ray volumes.
    n_sfc = wave_modes
    if nxray > 1
        n_sfc *= div(nxray, multiplication_factor)
    end
    if nyray > 1
        n_sfc *= div(nyray, multiplication_factor)
    end
    if nzray > 1
        n_sfc *= div(nzray, multiplication_factor)
    end

    # Allocate ray-volume arrays.
    nray = zeros(Int, nxx, nyy, nzz)
    rays = Rays(nray_wrk, nxx, nyy, nzz)
    merged_rays = MergedRays(2, nray_max)
    surface_indices = SurfaceIndices(n_sfc, nxx, nyy)
    increments = WKBIncrements(nray_wrk, nxx, nyy, nzz)
    integrals = WKBIntegrals(nxx, nyy, nzz)
    tendencies = WKBTendencies(nxx, nyy, nzz)

    spec_tend = TriadTendencies(nxx, nyy, nzz, kp_size, m_size)
    cgx_max = Ref(0.0)
    cgy_max = Ref(0.0)
    cgz_max = zeros(nxx, nyy, nzz)
    zb = zeros(nxx, nyy)
    diffusion = zeros(nxx, nyy, nzz)

    return WKB(
        nxray,
        nyray,
        nzray,
        nxray_wrk,
        nyray_wrk,
        nzray_wrk,
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        merged_rays,
        surface_indices,
        increments,
        integrals,
        tendencies,
        cgx_max,
        cgy_max,
        cgz_max,
        zb,
        diffusion,
        spec_tend
    )
end
