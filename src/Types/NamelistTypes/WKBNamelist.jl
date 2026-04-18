"""
```julia
WKBNamelist
```

Namelist for parameters used by MS-GWaM.

```julia
WKBNamelist(;
    nrx::Integer = 1,
    nry::Integer = 1,
    nrz::Integer = 1,
    nrk::Integer = 1,
    nrl::Integer = 1,
    nrm::Integer = 1,
    multiplication_factor::Integer = 4,
    dkr_factor::Real = 1.0E-1,
    dlr_factor::Real = 1.0E-1,
    dmr_factor::Real = 1.0E-1,
    branch::Integer = -1,
    merge_mode::Symbol = :ConstantWaveAction,
    filter_order::Integer = 2,
    smooth_tendencies::Bool = true,
    filter_type::Symbol = :ShapiroFilter,
    impact_altitude::Real = 0.0E+0,
    use_saturation::Bool = true,
    saturation_threshold::Real = 1.0E+0,
    wkb_mode::Symbol = :NoWKB,
    blocking::Bool = false,
    long_threshold::Real = 2.5E-1,
    drag_coefficient::Real = 1.0E+0,
    wave_modes::Integer = 1,
    initial_wave_field::Function = (alpha, x, y, z) ->
        (0.0, 0.0, 0.0, 0.0, 0.0),
    elastic_mode_selection::Bool = false,
    minimum_mode_count::Integer = wave_modes,
    maximum_mode_count::Integer = wave_modes,
    minimum_power_fraction::Real = 1.0E+0,
    maximum_power_fraction::Real = 1.0E+0,
)::WKBNamelist
```

Construct a `WKBNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `nrx::Int`: Number of ray-volumes launched per grid cell and wave mode in ``\\hat{x}``-direction.

  - `nry::Int`: Number of ray-volumes launched per grid cell and wave mode in ``\\hat{y}``-direction.

  - `nrz::Int`: Number of ray-volumes launched per grid cell and wave mode in ``\\hat{z}``-direction.

  - `nrk::Int`: Number of ray-volumes launched per grid cell and wave mode in ``k``-direction.

  - `nrl::Int`: Number of ray-volumes launched per grid cell and wave mode in ``l``-direction.

  - `nrm::Int`: Number of ray-volumes launched per grid cell and wave mode in ``m``-direction.

  - `multiplication_factor::Int`: Factor by which ray volumes are allowed to multiply in each dimension of physical space.

  - `dkr_factor::Float64`: Relative initial ray-volume extent in ``k``.

  - `dlr_factor::Float64`: Relative initial ray-volume extent in ``l``.

  - `dmr_factor::Float64`: Relative initial ray-volume extent in ``m``.

  - `branch::Int`: Frequency branch.

  - `merge_mode::Symbol`: Ray-volume merging strategy (conserved quantity).

  - `filter_order::Int`: Order of the smoothing applied to the mean-flow tendencies.

  - `smooth_tendencies::Bool`: Switch for smoothing the mean-flow tendencies.

  - `filter_type::Symbol`: Filter to use for the smoothing of the mean-flow tendencies.

  - `impact_altitude::Float64`: Minimum altitude for ray-tracing and mean-flow impact.

  - `use_saturation::Bool`: Switch for the saturation scheme.

  - `saturation_threshold::Float64`: Relative saturation threshold.

  - `wkb_mode::Symbol`: Approximations used by MS-GWaM.

  - `blocking::Bool`: Switch for parameterizing blocking in WKB-mountain-wave simulations.

  - `long_threshold::Float64`: Long-number threshold used by the blocked-layer scheme.

  - `drag_coefficient::Float64`: Dimensionless drag coefficient used by the blocked-layer scheme.

  - `wave_modes::Int`: Number of wave modes per grid cell.

  - `initial_wave_field::FunctionWrapper{NTuple{5, Float64}, Tuple{Int, Float64, Float64, Float64}}`: Function used to set the initial wavenumbers, intrinsic frequency and wave-action density of each wave mode.

  - `elastic_mode_selection::Bool`: Switch for elastic mode selection in ray-volume sources.

  - `minimum_mode_count::Int`: Minimum number of modes selected by the elastic-mode-selection algorithm.

  - `minimum_mode_count::Int`: Maximum number of modes selected by the elastic-mode-selection algorithm.

  - `minimum_power_fraction::Float64`: Minimum power fraction retained by the elastic-mode-selection algorithm.

  - `maximum_power_fraction::Float64`: Maximum power fraction retained by the elastic-mode-selection algorithm.

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been validated yet.
"""
struct WKBNamelist
    nrx::Int
    nry::Int
    nrz::Int
    nrk::Int
    nrl::Int
    nrm::Int
    multiplication_factor::Int
    dkr_factor::Float64
    dlr_factor::Float64
    dmr_factor::Float64
    branch::Int
    merge_mode::Symbol
    filter_order::Int
    smooth_tendencies::Bool
    filter_type::Symbol
    impact_altitude::Float64
    use_saturation::Bool
    saturation_threshold::Float64
    wkb_mode::Symbol
    blocking::Bool
    long_threshold::Float64
    drag_coefficient::Float64
    wave_modes::Int
    initial_wave_field::FunctionWrapper{
        NTuple{5, Float64},
        Tuple{Int, Float64, Float64, Float64},
    }
    elastic_mode_selection::Bool
    minimum_mode_count::Int
    maximum_mode_count::Int
    minimum_power_fraction::Float64
    maximum_power_fraction::Float64
end

function WKBNamelist(;
    nrx::Integer = 1,
    nry::Integer = 1,
    nrz::Integer = 1,
    nrk::Integer = 1,
    nrl::Integer = 1,
    nrm::Integer = 1,
    multiplication_factor::Integer = 4,
    dkr_factor::Real = 1.0E-1,
    dlr_factor::Real = 1.0E-1,
    dmr_factor::Real = 1.0E-1,
    branch::Integer = -1,
    merge_mode::Symbol = :ConstantWaveAction,
    filter_order::Integer = 2,
    smooth_tendencies::Bool = true,
    filter_type::Symbol = :ShapiroFilter,
    impact_altitude::Real = 0.0E+0,
    use_saturation::Bool = true,
    saturation_threshold::Real = 1.0E+0,
    wkb_mode::Symbol = :NoWKB,
    blocking::Bool = false,
    long_threshold::Real = 2.5E-1,
    drag_coefficient::Real = 1.0E+0,
    wave_modes::Integer = 1,
    initial_wave_field::Function = (alpha, x, y, z) ->
        (0.0, 0.0, 0.0, 0.0, 0.0),
    elastic_mode_selection::Bool = false,
    minimum_mode_count::Integer = wave_modes,
    maximum_mode_count::Integer = wave_modes,
    minimum_power_fraction::Real = 1.0E+0,
    maximum_power_fraction::Real = 1.0E+0,
)::WKBNamelist
    return WKBNamelist(
        Int(nrx),
        Int(nry),
        Int(nrz),
        Int(nrk),
        Int(nrl),
        Int(nrm),
        Int(multiplication_factor),
        Float64(dkr_factor),
        Float64(dlr_factor),
        Float64(dmr_factor),
        Int(branch),
        merge_mode,
        Int(filter_order),
        smooth_tendencies,
        filter_type,
        Float64(impact_altitude),
        use_saturation,
        Float64(saturation_threshold),
        wkb_mode,
        blocking,
        Float64(long_threshold),
        Float64(drag_coefficient),
        Int(wave_modes),
        initial_wave_field,
        elastic_mode_selection,
        Int(minimum_mode_count),
        Int(maximum_mode_count),
        Float64(minimum_power_fraction),
        Float64(maximum_power_fraction),
    )
end
