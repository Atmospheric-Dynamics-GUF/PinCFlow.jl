"""
```julia
WKBNamelist{
    A <: Integer,
    B <: AbstractFloat,
    C <: AbstractMergeMode,
    D <: Bool,
    E <: AbstractWKBFilter,
    F <: AbstractWKBMode,
    G <: Function,
}
```

Namelist for parameters used by MSGWaM.

```julia
WKBNamelist(;
    nrx::Integer = 1,
    nry::Integer = 1,
    nrz::Integer = 1,
    nrk::Integer = 1,
    nrl::Integer = 1,
    nrm::Integer = 1,
    multiplication_factor::Integer = 4,
    dkr_factor::AbstractFloat = 1.0E-1,
    dlr_factor::AbstractFloat = 1.0E-1,
    dmr_factor::AbstractFloat = 1.0E-1,
    branch::Integer = -1,
    merge_mode::AbstractMergeMode = ConstantWaveAction(),
    filter_order::Integer = 2,
    smooth_tendencies::Bool = true,
    filter_type::AbstractWKBFilter = Shapiro(),
    impact_altitude::AbstractFloat = 0.0E+0,
    use_saturation::Bool = true,
    saturation_threshold::AbstractFloat = 1.0E+0,
    wkb_mode::AbstractWKBMode = NoWKB(),
    blocking::Bool = false,
    long_threshold::AbstractFloat = 2.5E-1,
    drag_coefficient::AbstractFloat = 1.0E+0,
    wave_modes::Integer = 1,
    initial_wave_field::Function = (alpha, x, y, z) ->
        (0.0, 0.0, 0.0, 0.0, 0.0),
)::WKBNamelist
```

Construct a `WKBNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `nrx::A`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{x}``-direction.

  - `nry::A`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{y}``-direction.

  - `nrz::A`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{z}``-direction.

  - `nrk::A`: Number of ray-volumes launched per grid cell and wave mode in ``k``-direction.

  - `nrl::A`: Number of ray-volumes launched per grid cell and wave mode in ``l``-direction.

  - `nrm::A`: Number of ray-volumes launched per grid cell and wave mode in ``m``-direction.

  - `multiplication_factor::A`: Factor by which ray volumes are allowed to multiply in each dimension of physical space.

  - `dkr_factor::B`: Relative initial ray-volume extent in ``k``.

  - `dlr_factor::B`: Relative initial ray-volume extent in ``l``.

  - `dmr_factor::B`: Relative initial ray-volume extent in ``m``.

  - `branch::A`: Frequency branch.

  - `merge_mode::C`: Ray-volume merging strategy (conserved quantity).

  - `filter_order::A`: Order of the smoothing applied to the mean-flow tendencies.

  - `smooth_tendencies::D`: Switch for smoothing the mean-flow tendencies.

  - `filter_type::E`: Filter to use for the smoothing of the mean-flow tendencies.

  - `impact_altitude::B`: Minimum altitude for ray-tracing and mean-flow impact.

  - `use_saturation::D`: Switch for the saturation scheme.

  - `saturation_threshold::B`: Relative saturation threshold.

  - `wkb_mode::F`: Approximations used by MSGWaM.

  - `blocking::D`: Switch for parameterizing blocking in WKB-mountain-wave simulations.

  - `long_threshold::B`: Long-number threshold used by the blocked-layer scheme.

  - `drag_coefficient::B`: Dimensionless (relative) drag coefficient used by the blocked layer scheme.

  - `wave_modes::A`: Number of wave modes per grid cell.

  - `initial_wave_field::G`: Function used to set the initial wavenumbers, intrinsic frequency and wave-action density of each wave mode.
"""
struct WKBNamelist{
    A <: Integer,
    B <: AbstractFloat,
    C <: AbstractMergeMode,
    D <: Bool,
    E <: AbstractWKBFilter,
    F <: AbstractWKBMode,
    G <: Function,
}
    nrx::A
    nry::A
    nrz::A
    nrk::A
    nrl::A
    nrm::A
    multiplication_factor::A
    dkr_factor::B
    dlr_factor::B
    dmr_factor::B
    branch::A
    merge_mode::C
    filter_order::A
    smooth_tendencies::D
    filter_type::E
    impact_altitude::B
    use_saturation::D
    saturation_threshold::B
    wkb_mode::F
    blocking::D
    long_threshold::B
    drag_coefficient::B
    wave_modes::A
    initial_wave_field::G
end

function WKBNamelist(;
    nrx::Integer = 1,
    nry::Integer = 1,
    nrz::Integer = 1,
    nrk::Integer = 1,
    nrl::Integer = 1,
    nrm::Integer = 1,
    multiplication_factor::Integer = 4,
    dkr_factor::AbstractFloat = 1.0E-1,
    dlr_factor::AbstractFloat = 1.0E-1,
    dmr_factor::AbstractFloat = 1.0E-1,
    branch::Integer = -1,
    merge_mode::AbstractMergeMode = ConstantWaveAction(),
    filter_order::Integer = 2,
    smooth_tendencies::Bool = true,
    filter_type::AbstractWKBFilter = Shapiro(),
    impact_altitude::AbstractFloat = 0.0E+0,
    use_saturation::Bool = true,
    saturation_threshold::AbstractFloat = 1.0E+0,
    wkb_mode::AbstractWKBMode = NoWKB(),
    blocking::Bool = false,
    long_threshold::AbstractFloat = 2.5E-1,
    drag_coefficient::AbstractFloat = 1.0E+0,
    wave_modes::Integer = 1,
    initial_wave_field::Function = (alpha, x, y, z) ->
        (0.0, 0.0, 0.0, 0.0, 0.0),
)::WKBNamelist
    return WKBNamelist(
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        multiplication_factor,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        branch,
        merge_mode,
        filter_order,
        smooth_tendencies,
        filter_type,
        impact_altitude,
        use_saturation,
        saturation_threshold,
        wkb_mode,
        blocking,
        long_threshold,
        drag_coefficient,
        wave_modes,
        initial_wave_field,
    )
end
