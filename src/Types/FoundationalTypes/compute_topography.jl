"""
```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    test_case::AbstractWKBTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
```

Compute and return the topography for WKB test cases.

The arrays in the returned tuple represent (in order) the resolved topography, the amplitudes of the unresolved topography, the corresponding zonal wavenumbers and the corresponding meridional wavenumbers.

```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    test_case::AbstractTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
```

Compute and return the topography for non-WKB test cases.

The arrays representing the unresolved spectrum are set to have the size `(0, 0, 0)`. The topography is represented by the first array in the returned tuple.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `x`: ``\\widehat{x}``-coordinate grid points.

  - `y`: ``\\widehat{y}``-coordinate grid points.

  - `test_case`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_topography end

function compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    test_case::AbstractWKBTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
    (; wave_modes) = namelists.wkb
    (; spectral_modes, resolved_topography, unresolved_topography) =
        namelists.grid
    (; nxx, nyy, io, jo, i0, i1, j0, j1) = domain
    (; lref) = constants

    if wave_modes < 1 || wave_modes < spectral_modes
        error("Error in compute_topography: wave_modes is too small!")
    end

    hb = zeros(nxx, nyy)
    hw = zeros(wave_modes, nxx, nyy)
    kh = zeros(wave_modes, nxx, nyy)
    lh = zeros(wave_modes, nxx, nyy)

    @ivy for j in j0:j1, i in i0:i1
        hbdim =
            resolved_topography(namelists, x[io + i] * lref, y[jo + j] * lref)
        hb[i, j] = hbdim / lref
        for alpha in 1:spectral_modes
            (khdim, lhdim, hwdim) = unresolved_topography(
                namelists,
                alpha,
                x[io + i] * lref,
                y[jo + j] * lref,
            )
            kh[alpha, i, j] = khdim * lref
            lh[alpha, i, j] = lhdim * lref
            hw[alpha, i, j] = hwdim / lref
        end
    end

    set_zonal_boundaries_of_field!(hb, namelists, domain)
    set_meridional_boundaries_of_field!(hb, namelists, domain)

    return (hb, hw, kh, lh)
end

function compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    test_case::AbstractTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
    (; resolved_topography) = namelists.grid
    (; nxx, nyy, io, jo, i0, i1, j0, j1) = domain
    (; lref) = constants

    hb = zeros(nxx, nyy)
    hw = zeros(0, 0, 0)
    kh = zeros(0, 0, 0)
    lh = zeros(0, 0, 0)

    @ivy for j in j0:j1, i in i0:i1
        hbdim =
            resolved_topography(namelists, x[io + i] * lref, y[jo + j] * lref)
        hb[i, j] = hbdim / lref
    end

    set_zonal_boundaries_of_field!(hb, namelists, domain)
    set_meridional_boundaries_of_field!(hb, namelists, domain)

    return (hb, hw, kh, lh)
end
