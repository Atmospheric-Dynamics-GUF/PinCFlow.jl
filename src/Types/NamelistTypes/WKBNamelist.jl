"""
```julia
WKBNamelist{
    A <: Integer,
    B <: Integer,
    C <: Integer,
    D <: Integer,
    E <: Integer,
    F <: Integer,
    G <: Integer,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: Integer,
    L <: AbstractMergeMode,
    M <: Integer,
    N <: Bool,
    O <: AbstractWKBFilter,
    P <: AbstractFloat,
    Q <: Bool,
    R <: AbstractFloat,
    S <: AbstractWKBMode,
    T <: Bool,
    U <: AbstractFloat,
    V <: AbstractFloat,
    W <: Integer,
    X <: Function,
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

  - `nry::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{y}``-direction.

  - `nrz::C`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{z}``-direction.

  - `nrk::D`: Number of ray-volumes launched per grid cell and wave mode in ``k``-direction.

  - `nrl::E`: Number of ray-volumes launched per grid cell and wave mode in ``l``-direction.

  - `nrm::F`: Number of ray-volumes launched per grid cell and wave mode in ``m``-direction.

  - `multiplication_factor::G`: Factor by which ray volumes are allowed to multiply in each dimension of physical space.

  - `dkr_factor::H`: Relative initial ray-volume extent in ``k``.

  - `dlr_factor::I`: Relative initial ray-volume extent in ``l``.

  - `dmr_factor::J`: Relative initial ray-volume extent in ``m``.

  - `branch::K`: Frequency branch.

  - `merge_mode::L`: Ray-volume merging strategy (conserved quantity).

  - `filter_order::M`: Order of the smoothing applied to the mean-flow tendencies.

  - `smooth_tendencies::N`: Switch for smoothing the mean-flow tendencies.

  - `filter_type::O`: Filter to use for the smoothing of the mean-flow tendencies.

  - `impact_altitude::P`: Minimum altitude for ray-tracing and mean-flow impact.

  - `use_saturation::Q`: Switch for the saturation scheme.

  - `saturation_threshold::R`: Relative saturation threshold.

  - `wkb_mode::S`: Approximations used by MSGWaM.

  - `blocking::T`: Switch for parameterizing blocking in WKB-mountain-wave simulations.

  - `long_threshold::U`: Long-number threshold used by the blocked-layer scheme.

  - `drag_coefficient::V`: Dimensionless (relative) drag coefficient used by the blocked layer scheme.

  - `wave_modes::W`: Number of wave modes per grid cell.

  - `initial_wave_field::X`: Function used to set the initial wavenumbers, intrinsic frequency and wave-action density of each wave mode.
"""
struct WKBNamelist{
    A <: Integer,
    B <: Integer,
    C <: Integer,
    D <: Integer,
    E <: Integer,
    F <: Integer,
    G <: Integer,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: Integer,
    L <: AbstractMergeMode,
    M <: Integer,
    N <: Bool,
    O <: AbstractWKBFilter,
    P <: AbstractFloat,
    Q <: Bool,
    R <: AbstractFloat,
    S <: AbstractWKBMode,
    T <: Bool,
    U <: AbstractFloat,
    V <: AbstractFloat,
    W <: Integer,
    X <: Function,
}
    nrx::A
    nry::B
    nrz::C
    nrk::D
    nrl::E
    nrm::F
    multiplication_factor::G
    dkr_factor::H
    dlr_factor::I
    dmr_factor::J
    branch::K
    merge_mode::L
    filter_order::M
    smooth_tendencies::N
    filter_type::O
    impact_altitude::P
    use_saturation::Q
    saturation_threshold::R
    wkb_mode::S
    blocking::T
    long_threshold::U
    drag_coefficient::V
    wave_modes::W
    initial_wave_field::X
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
