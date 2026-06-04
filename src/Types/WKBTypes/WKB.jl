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
    J <: AbstractMatrix{<:AbstractFloat},
    K <: AbstractArray{<:AbstractFloat, 3},
}
```

Main container for WKB ray-tracing data and parameters.

```julia
WKB(namelists::Namelists, domain::Domain)::WKB
```

Construct a `WKB` instance by dispatching to a test-case-specific method.

```julia
WKB(namelists::Namelists, domain::Domain, wkb_mode::Val{:NoWKB})::WKB
```

Construct a `WKB` instance with zero-size arrays for non-WKB configurations.

```julia
WKB(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::WKB
```

Construct a `WKB` instance.

This method primarily determines the size of the spectral dimension of ray-volume arrays and initializes them and related arrays (with zeros) accordingly. The proper initialization with nonzero wave action is performed by [`PinCFlow.MSGWaM.RayUpdate.initialize_rays!`](@ref).

# Fields

  - `bins::A`: Maximum ray-volume count allowed per grid cell before merging is triggered (equal to `n_sfc` in steady-state mode and `k_bins * l_bins * m_bins` otherwise).

  - `nray_wrk::A`: Size of the spectral dimension of ray-volume arrays.

  - `n_sfc::A`: Number of orographic wave modes.

  - `nray::B`: Ray-volume count in each grid cell.

  - `rays::C`: Prognostic ray-volume properties.

  - `merged_rays::D`: Container used for creating merged ray volumes.

  - `surface_indices::E`: Indices that connect orographic wave modes to ray volumes.

  - `increments::F`: WKBIncrements of the prognostic ray-volume properties.

  - `integrals::G`: Integrals of ray-volume properties.

  - `tendencies::H`: Gravity-wave drag and heating fields.

  - `cgx_max::I`: Maximum zonal group velocity.

  - `cgy_max::I`: Maximum meridional group velocity.

  - `cgz_max::I`: Maximum vertical group velocity.

  - `zb::J`: Upper edge of the blocked layer.

  - `diffusion::K`: Diffusion induced by wave breaking.

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
    J <: AbstractMatrix{<:AbstractFloat},
    K <: AbstractArray{<:AbstractFloat, 3},
}
    bins::A
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
    cgz_max::I
    zb::J
    diffusion::K
end

function WKB(namelists::Namelists, domain::Domain)::WKB
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return WKB(namelists, domain, Val(wkb_mode))
end

function WKB(namelists::Namelists, domain::Domain, wkb_mode::Val{:NoWKB})::WKB
    return WKB(
        [0 for i in 1:3]...,
        zeros(Int, 0, 0, 0),
        Rays(0, 0, 0, 0),
        MergedRays(0, 0),
        SurfaceIndices(0, 0, 0),
        WKBIncrements(0, 0, 0, 0),
        WKBIntegrals(0, 0, 0),
        WKBTendencies(0, 0, 0),
        [Ref(0.0) for i in 1:3]...,
        zeros(0, 0),
        zeros(0, 0, 0),
    )
end

function WKB(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::WKB
    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        k_bins,
        l_bins,
        m_bins,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
    ) = namelists.wkb
    (; x_size, y_size, z_size) = namelists.domain
    (; nxx, nyy, nzz) = domain

    # Check if the spectral-extent factors are set correctly.
    if x_size > 1 && dkr_factor == 0.0
        error("Error in WKB: x_size > 1 && dkr_factor == 0!")
    end
    if y_size > 1 && dlr_factor == 0.0
        error("Error in WKB: y_size > 1 && dlr_factor == 0!")
    end
    if z_size == 1 || dmr_factor == 0.0
        error("Error in WKB: z_size == 1 || dmr_factor == 0!")
    end

    # Check if the numbers of bins are set correctly.
    if k_bins % 2 == 0
        error("k_bins must be an odd number!")
    end
    if l_bins % 2 == 0
        error("l_bins must be an odd number!")
    end
    if m_bins % 2 == 0
        error("m_bins must be an odd number!")
    end

    # Set the number of surface ray volumes.
    n_sfc = nrx * nrk * nry * nrl * nrz * nrm * wave_modes

    # Set the total number of bins and work size of the ray-volume array.
    if wkb_mode == Val(:SteadyState)
        bins = nray_wrk = n_sfc
    else
        # Set the total number of bins.
        bins = k_bins * l_bins * m_bins

        # Check if the number of bins is large enough.
        if bins < n_sfc
            error(
                "k_bins * l_bins * m_bins must not be smaller than nrx * nrk * nry * nrl * nrz * nrm * wave_modes",
            )
        end

        # Determine the work size of the ray-volume array.
        nray_wrk = 2 * bins
        y_size > 1 && (nray_wrk *= 2)
        x_size > 1 && (nray_wrk *= 2)
    end

    # Allocate ray-volume arrays.
    nray = zeros(Int, nxx, nyy, nzz)
    rays = Rays(nray_wrk, nxx, nyy, nzz)
    merged_rays = MergedRays(2, bins)
    surface_indices = SurfaceIndices(n_sfc, nxx, nyy)
    increments = WKBIncrements(nray_wrk, nxx, nyy, nzz)
    integrals = WKBIntegrals(nxx, nyy, nzz)
    tendencies = WKBTendencies(nxx, nyy, nzz)
    cgx_max = Ref(0.0)
    cgy_max = Ref(0.0)
    cgz_max = Ref(0.0)
    zb = zeros(nxx, nyy)
    diffusion = zeros(nxx, nyy, nzz)

    return WKB(
        bins,
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
    )
end
