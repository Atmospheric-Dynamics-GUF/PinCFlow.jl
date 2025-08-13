"""
```julia
WKB{
    A <: Integer,
    B <: AbstractArray{<:Integer, 3},
    C <: Rays,
    D <: SurfaceIndices,
    E <: Increments,
    F <: GWIntegrals,
    G <: GWTendencies,
    H <: Ref{<:AbstractFloat},
    I <: AbstractArray{<:AbstractFloat, 3},
    J <: AbstractMatrix{<:AbstractFloat},
}
```

Main container for WKB ray-tracing data and parameters.

```julia
WKB(namelists::Namelists, constants::Constants, domain::Domain, grid::Grid)
```

Construct a `WKB` instance by dispatching to a test-case-specific method.

```julia
WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    testcase::AbstractTestCase,
)
```

Construct a `WKB` instance with zero-size arrays for non-WKB test cases.

```julia
WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    testcase::AbstractWKBTestCase,
)
```

Construct a `WKB` instance.

This method primarily determines the size of the spectral dimension of ray-volume arrays and initializes them and related arrays (with zeros) accordingly. The proper initialization with nonzero wave action is performed by [`PinCFlow.MSGWaM.RayUpdate.initialize_rays!`](@ref).

# Fields

  - `nxray::A`: Number of ray volumes allowed in ``\\widehat{x}``, per grid cell and wave mode (`nray_fac * nrxl * nrk_init`, taken from `namelists.wkb`).

  - `nyray::A`: Number of ray volumes allowed in ``\\widehat{y}``, per grid cell and wave mode (`nray_fac * nryl * nrl_init`, taken from `namelists.wkb`).

  - `nzray::A`: Number of ray volumes allowed in ``\\widehat{z}``, per grid cell and wave mode (`nray_fac * nrzl * nrm_init`, taken from `namelists.wkb`).

  - `nxray_wrk::A`: `2 * nxray`.

  - `nyray_wrk::A`: `2 * nyray`.

  - `nzray_wrk::A`: `2 * nzray`.

  - `nray_max::A`: Maximum ray-volume count allowed per grid-cell before merging is triggered (`nxray * nyray * nzray * namelists.wkb.nwm`).

  - `nray_wrk::A`: Size of the spectral dimension of ray-volume arrays (`nxray_wrk * nyray_wrk * nzray_wrk`).

  - `n_sfc::A`: Number of orographic wave modes.

  - `nray::B`: Ray-volume count in each grid cell.

  - `rays::C`: Prognostic ray-volume properties.

  - `surface_indices::D`: Indices that connect orographic wave modes to ray volumes.

  - `increments::E`: Increments of the prognostic ray-volume properties.

  - `integrals::F`: Integrals of ray-volume properties.

  - `tendencies::G`: Gravity-wave drag and heating fields.

  - `cgx_max::H`: Maximum zonal group velocities.

  - `cgy_max::H`: Maximum meridional group velocities.

  - `cgz_max::I`: Maximum vertical group velocities.

  - `zb::J`: Upper edge of the blocked layer.

  - `diffusion::I`: Diffusion induced by wave breaking.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `testcase`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.Types.WKBTypes.Rays`](@ref)

  - [`PinCFlow.Types.WKBTypes.SurfaceIndices`](@ref)

  - [`PinCFlow.Types.WKBTypes.Increments`](@ref)

  - [`PinCFlow.Types.WKBTypes.GWIntegrals`](@ref)

  - [`PinCFlow.Types.WKBTypes.GWTendencies`](@ref)
"""
struct WKB{
    A <: Integer,
    B <: AbstractArray{<:Integer, 3},
    C <: Rays,
    D <: SurfaceIndices,
    E <: Increments,
    F <: GWIntegrals,
    G <: GWTendencies,
    H <: Ref{<:AbstractFloat},
    I <: AbstractArray{<:AbstractFloat, 3},
    J <: AbstractMatrix{<:AbstractFloat},
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
    surface_indices::D
    increments::E
    integrals::F
    tendencies::G
    cgx_max::H
    cgy_max::H
    cgz_max::I
    zb::J
    diffusion::I
end

function WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)
    (; testcase) = namelists.setting
    return WKB(namelists, constants, domain, grid, testcase)
end

function WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    testcase::AbstractTestCase,
)
    return WKB(
        [0 for i in 1:9]...,
        zeros(Int, 0, 0, 0),
        Rays(0, 0, 0, 0),
        SurfaceIndices(0, 0, 0),
        Increments(0, 0, 0, 0),
        GWIntegrals(0, 0, 0),
        GWTendencies(0, 0, 0),
        [Ref(0.0) for i in 1:2]...,
        zeros(0, 0, 0),
        zeros(0, 0),
        zeros(0, 0, 0),
    )
end

function WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    testcase::AbstractWKBTestCase,
)
    (;
        xrmin_dim,
        xrmax_dim,
        yrmin_dim,
        yrmax_dim,
        zrmin_dim,
        zrmax_dim,
        nrxl,
        nryl,
        nrzl,
        nrk_init,
        nrl_init,
        nrm_init,
        nray_fac,
        nwm,
        fac_dk_init,
        fac_dl_init,
        fac_dm_init,
    ) = namelists.wkb
    (; lref) = constants
    (; sizex, sizey, sizez) = namelists.domain
    (; nxx, nyy, nzz) = domain
    (; lx, ly, lz) = grid

    # Non-dimensionalize boundaries for ray-volume propagation.
    xrmin = xrmin_dim / lref
    xrmax = xrmax_dim / lref
    yrmin = yrmin_dim / lref
    yrmax = yrmax_dim / lref
    zrmin = zrmin_dim / lref
    zrmax = zrmax_dim / lref

    # Check if the boundaries for ray-volume propagation are within the domain.
    if xrmin < lx[1] || xrmax > lx[2]
        error("Error in WKB: xrmin too small or xrmax too large!")
    end
    if yrmin < ly[1] || yrmax > ly[2]
        error("Error in WKB: yrmin too small or yrmax too large!")
    end
    if zrmin < lz[1] || zrmax > lz[2]
        error("Error in WKB: zrmin too small or zrmax too large!")
    end

    # Check if spectral-extent factors are set correctly.
    if sizex > 1 && fac_dk_init == 0.0
        error("Error in WKB: sizex > 1 && fac_dk_init == 0!")
    end
    if sizey > 1 && fac_dl_init == 0.0
        error("Error in WKB: sizey > 1 && fac_dl_init == 0!")
    end
    if sizez == 1 || fac_dm_init == 0.0
        error("Error in WKB: sizez == 1 || fac_dm_init == 0!")
    end

    # Set zonal ray-volume count.
    if sizex == 1
        nxray = 1
    else
        nxray = nray_fac * nrxl * nrk_init
    end

    # Set meridional ray-volume count.
    if sizey == 1
        nyray = 1
    else
        nyray = nray_fac * nryl * nrl_init
    end

    # Set vertical ray-volume count.
    nzray = nray_fac * nrzl * nrm_init

    # Set maximum ray-volume count.
    nray_max = nxray * nyray * nzray * nwm

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
    nray_wrk = nxray_wrk * nyray_wrk * nzray_wrk

    # Set number of surface ray volumes.
    n_sfc = nwm
    if nxray > 1
        n_sfc *= div(nxray, nray_fac)
    end
    if nyray > 1
        n_sfc *= div(nyray, nray_fac)
    end
    if nzray > 1
        n_sfc *= div(nzray, nray_fac)
    end

    # Initialize ray-volume arrays.
    nray = zeros(Int, nxx, nyy, nzz)
    rays = Rays(nray_wrk, nxx, nyy, nzz)
    surface_indices = SurfaceIndices(n_sfc, nxx, nyy)
    increments = Increments(nray_wrk, nxx, nyy, nzz)
    integrals = GWIntegrals(nxx, nyy, nzz)
    tendencies = GWTendencies(nxx, nyy, nzz)
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
        surface_indices,
        increments,
        integrals,
        tendencies,
        cgx_max,
        cgy_max,
        cgz_max,
        zb,
        diffusion,
    )
end
