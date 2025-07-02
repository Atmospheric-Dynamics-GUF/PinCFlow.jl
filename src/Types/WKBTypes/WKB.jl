"""
    WKB

Main container for WKB ray tracing data and parameters.

# Fields

  - `nxray`, `nyray`, `nzray`: Number of ray volumes in each direction
  - `nxray_wrk`, `nyray_wrk`, `nzray_wrk`: Working array dimensions (spectral)
  - `nray_max`, `nray_wrk`: Maximum and working ray counts
  - `n_sfc`: Number of surface ray volumes
  - `nray`: 3D array of ray counts per grid cell
  - `rays`: Ray position, wavenumber, and density data
  - `surface_indices`: Indices for surface ray launching
  - `increments`: Ray propagation increments
  - `integrals`: Gravity wave momentum and energy integrals
  - `tendencies`: Gravity wave drag and heating tendencies
  - `cgx_max`, `cgy_max`, `cgz_max`: Maximum group velocities
  - `zb`: Bottom boundary height
  - `diffusion`: Diffusion coefficients
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
