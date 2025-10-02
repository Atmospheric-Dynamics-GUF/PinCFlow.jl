"""
```julia
WKBNamelist{
    A <: AbstractFloat,
    B <: Integer,
    C <: AbstractMergeMode,
    D <: Bool,
    E <: AbstractWKBFilter,
    F <: AbstractWKBMode,
}
```

Namelist for parameters used by MSGWaM.

```julia
WKBNamelist(;
    xrmin::AbstractFloat = 0.0E+0,
    xrmax::AbstractFloat = 1.0E+3,
    yrmin::AbstractFloat = 0.0E+0,
    yrmax::AbstractFloat = 1.0E+3,
    zrmin::AbstractFloat = 0.0E+0,
    zrmax::AbstractFloat = 1.0E+3,
    nrx::Integer = 1,
    nry::Integer = 1,
    nrz::Integer = 1,
    nrk::Integer = 1,
    nrl::Integer = 1,
    nrm::Integer = 1,
    multiplication_factor::Integer = 4,
    fdk::AbstractFloat = 1.0E-1,
    fdl::AbstractFloat = 1.0E-1,
    fdm::AbstractFloat = 1.0E-1,
    branch::Integer = -1,
    merge_mode::AbstractMergeMode = ConstantWaveAction(),
    filter_order::Integer = 2,
    smooth_tendencies::Bool = true,
    filter_type::AbstractWKBFilter = Shapiro(),
    impact_altitude::AbstractFloat = 0.0E+0,
    use_saturation::Bool = true,
    saturation_threshold::AbstractFloat = 1.0E+0,
    wkb_mode::AbstractWKBMode = MultiColumn(),
    blocking::Bool = false,
    long_threshold::AbstractFloat = 2.5E-1,
    drag_coefficient::AbstractFloat = 1.0E+0,
    nalpha::Integer = 1,
)::WKBNamelist
```

Construct a `WKBNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `xrmin::A`: Lower bound in ``x`` of the region where ray volumes are initialized.

  - `xrmax::A`: Upper bound in ``x`` of the region where ray volumes are initialized.

  - `yrmin::A`: Lower bound in ``y`` of the region where ray volumes are initialized.

  - `yrmax::A`: Upper bound in ``y`` of the region where ray volumes are initialized.

  - `zrmin::A`: Lower bound in ``z`` of the region where ray volumes are initialized.

  - `zrmax::A`: Upper bound in ``z`` of the region where ray volumes are initialized.

  - `nrx::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{x}``-direction.

  - `nry::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{y}``-direction.

  - `nrz::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{z}``-direction.

  - `nrk::B`: Number of ray-volumes launched per grid cell and wave mode in ``k``-direction.

  - `nrl::B`: Number of ray-volumes launched per grid cell and wave mode in ``l``-direction.

  - `nrm::B`: Number of ray-volumes launched per grid cell and wave mode in ``m``-direction.

  - `multiplication_factor::B`: Factor by which ray volumes are allowed to multiply in each dimension of physical space.

  - `fdk::A`: Relative initial ray-volume extent in ``k``.

  - `fdl::A`: Relative initial ray-volume extent in ``l``.

  - `fdm::A`: Relative initial ray-volume extent in ``m``.

  - `branch::B`: Frequency branch.

  - `merge_mode::C`: Ray-volume merging strategy (conserved quantity).

  - `filter_order::B`: Order of the smoothing applied to the mean-flow tendencies.

  - `smooth_tendencies::D`: Switch for smoothing the mean-flow tendencies.

  - `filter_type::E`: Filter to use for the smoothing of the mean-flow tendencies.

  - `impact_altitude::A`: Minimum altitude for ray-tracing and mean-flow impact.

  - `use_saturation::D`: Switch for the saturation scheme.

  - `saturation_threshold::A`: Relative saturation threshold.

  - `wkb_mode::F`: Approximations used by MSGWaM.

  - `blocking::D`: Switch for parameterizing blocking in WKB-mountain-wave simulations.

  - `long_threshold::A`: Long-number threshold used by the blocked-layer scheme.

  - `drag_coefficient::A`: Dimensionless (relative) drag coefficient used by the blocked layer scheme.

  - `nalpha::B`: Number of wave modes per grid cell.
"""
struct WKBNamelist{
    A <: AbstractFloat,
    B <: Integer,
    C <: AbstractMergeMode,
    D <: Bool,
    E <: AbstractWKBFilter,
    F <: AbstractWKBMode,
}
    xrmin::A
    xrmax::A
    yrmin::A
    yrmax::A
    zrmin::A
    zrmax::A
    nrx::B
    nry::B
    nrz::B
    nrk::B
    nrl::B
    nrm::B
    multiplication_factor::B
    fdk::A
    fdl::A
    fdm::A
    branch::B
    merge_mode::C
    filter_order::B
    smooth_tendencies::D
    filter_type::E
    impact_altitude::A
    use_saturation::D
    saturation_threshold::A
    wkb_mode::F
    blocking::D
    long_threshold::A
    drag_coefficient::A
    nalpha::B
end

function WKBNamelist(;
    xrmin::AbstractFloat = 0.0E+0,
    xrmax::AbstractFloat = 1.0E+3,
    yrmin::AbstractFloat = 0.0E+0,
    yrmax::AbstractFloat = 1.0E+3,
    zrmin::AbstractFloat = 0.0E+0,
    zrmax::AbstractFloat = 1.0E+3,
    nrx::Integer = 1,
    nry::Integer = 1,
    nrz::Integer = 1,
    nrk::Integer = 1,
    nrl::Integer = 1,
    nrm::Integer = 1,
    multiplication_factor::Integer = 4,
    fdk::AbstractFloat = 1.0E-1,
    fdl::AbstractFloat = 1.0E-1,
    fdm::AbstractFloat = 1.0E-1,
    branch::Integer = -1,
    merge_mode::AbstractMergeMode = ConstantWaveAction(),
    filter_order::Integer = 2,
    smooth_tendencies::Bool = true,
    filter_type::AbstractWKBFilter = Shapiro(),
    impact_altitude::AbstractFloat = 0.0E+0,
    use_saturation::Bool = true,
    saturation_threshold::AbstractFloat = 1.0E+0,
    wkb_mode::AbstractWKBMode = MultiColumn(),
    blocking::Bool = false,
    long_threshold::AbstractFloat = 2.5E-1,
    drag_coefficient::AbstractFloat = 1.0E+0,
    nalpha::Integer = 1,
)::WKBNamelist
    return WKBNamelist(
        xrmin,
        xrmax,
        yrmin,
        yrmax,
        zrmin,
        zrmax,
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        multiplication_factor,
        fdk,
        fdl,
        fdm,
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
        nalpha,
    )
end
