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

Namelist for the configuration of MS-GWaM (see constructor for parameter descriptions).
"""
struct WKBNamelist{
    A <: AbstractFloat,
    B <: Integer,
    C <: AbstractMergeMode,
    D <: Bool,
    E <: AbstractWKBFilter,
    F <: AbstractWKBMode,
}
    xrmin_dim::A
    xrmax_dim::A
    yrmin_dim::A
    yrmax_dim::A
    zrmin_dim::A
    zrmax_dim::A
    nrxl::B
    nryl::B
    nrzl::B
    nrk_init::B
    nrl_init::B
    nrm_init::B
    nray_fac::B
    fac_dk_init::A
    fac_dl_init::A
    fac_dm_init::A
    branchr::B
    merge_mode::C
    nsmth_wkb::B
    lsmth_wkb::D
    sm_filter::E
    zmin_wkb_dim::A
    lsaturation::D
    alpha_sat::A
    wkb_mode::F
    blocking::D
    long_threshold::A
    drag_coefficient::A
    nwm::B
end

"""
```julia
WKBNamelist(;
    xrmin_dim = 0.0E+0,
    xrmax_dim = 1.0E+3,
    yrmin_dim = 0.0E+0,
    yrmax_dim = 1.0E+3,
    zrmin_dim = 0.0E+0,
    zrmax_dim = 1.0E+3,
    nrxl = 1,
    nryl = 1,
    nrzl = 1,
    nrk_init = 1,
    nrl_init = 1,
    nrm_init = 1,
    nray_fac = 4,
    fac_dk_init = 1.0E-1,
    fac_dl_init = 1.0E-1,
    fac_dm_init = 1.0E-1,
    branchr = -1,
    merge_mode = ConstantWaveAction(),
    nsmth_wkb = 2,
    lsmth_wkb = true,
    sm_filter = Shapiro(),
    zmin_wkb_dim = 0.0,
    lsaturation = true,
    alpha_sat = 1.0E+0,
    wkb_mode = MultiColumn(),
    blocking = false,
    long_threshold = 2.5E-1,
    drag_coefficient = 1.0E+0,
    nwm = 1,
)
```

Construct a WKBNamelist instance, which holds parameters for the WKB ray tracing algorithm.

# Arguments:

  - `xrmin_dim`: Minimum x-coordinate of ray tracing domain in dimensional units.
  - `xrmax_dim`: Maximum x-coordinate of ray tracing domain in dimensional units.
  - `yrmin_dim`: Minimum y-coordinate of ray tracing domain in dimensional units.
  - `yrmax_dim`: Maximum y-coordinate of ray tracing domain in dimensional units.
  - `zrmin_dim`: Minimum z-coordinate of ray tracing domain in dimensional units.
  - `zrmax_dim`: Maximum z-coordinate of ray tracing domain in dimensional units.
  - `nrxl`: Number of ray launch points in x-direction.
  - `nryl`: Number of ray launch points in y-direction.
  - `nrzl`: Number of ray launch points in z-direction.
  - `nrk_init`: Initial number of wavenumber points in x-direction.
  - `nrl_init`: Initial number of wavenumber points in y-direction.
  - `nrm_init`: Initial number of wavenumber points in z-direction.
  - `nray_fac`: Ray multiplication factor.
  - `fac_dk_init`: Initial wavenumber spacing factor in x-direction.
  - `fac_dl_init`: Initial wavenumber spacing factor in y-direction.
  - `fac_dm_init`: Initial wavenumber spacing factor in z-direction.
  - `branchr`: Branch selection for ray tracing.
  - `merge_mode`: Algorithm to use when merging ray beams.
  - `nsmth_wkb`: Number of smoothing passes for WKB solution.
  - `lsmth_wkb`: Whether to apply smoothing to WKB solution.
  - `sm_filter`: Type of smoothing filter to apply to WKB solution.
  - `zmin_wkb_dim`: Minimum height for WKB calculation in dimensional units.
  - `lsaturation`: Whether to apply wave saturation.
  - `alpha_sat`: Saturation coefficient.
  - `wkb_mode`: WKB calculation mode.
  - `blocking`: Whether to apply wave blocking.
  - `long_threshold`: Threshold for long wavelength approximation.
  - `drag_coefficient`: Coefficient for wave-induced drag.
  - `nwm`: Number of wave modes.
  - `launch_algorithm`: Algorithm used for wave launching.

# Returns

  - `::WKBNamelist`: `WKBNamelist` instance.
"""
function WKBNamelist(;
    xrmin_dim = 0.0E+0,
    xrmax_dim = 1.0E+3,
    yrmin_dim = 0.0E+0,
    yrmax_dim = 1.0E+3,
    zrmin_dim = 0.0E+0,
    zrmax_dim = 1.0E+3,
    nrxl = 1,
    nryl = 1,
    nrzl = 1,
    nrk_init = 1,
    nrl_init = 1,
    nrm_init = 1,
    nray_fac = 4,
    fac_dk_init = 1.0E-1,
    fac_dl_init = 1.0E-1,
    fac_dm_init = 1.0E-1,
    branchr = -1,
    merge_mode = ConstantWaveAction(),
    nsmth_wkb = 2,
    lsmth_wkb = true,
    sm_filter = Shapiro(),
    zmin_wkb_dim = 0.0,
    lsaturation = true,
    alpha_sat = 1.0E+0,
    wkb_mode = MultiColumn(),
    blocking = false,
    long_threshold = 2.5E-1,
    drag_coefficient = 1.0E+0,
    nwm = 1,
)
    return WKBNamelist(
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
        fac_dk_init,
        fac_dl_init,
        fac_dm_init,
        branchr,
        merge_mode,
        nsmth_wkb,
        lsmth_wkb,
        sm_filter,
        zmin_wkb_dim,
        lsaturation,
        alpha_sat,
        wkb_mode,
        blocking,
        long_threshold,
        drag_coefficient,
        nwm,
    )
end
