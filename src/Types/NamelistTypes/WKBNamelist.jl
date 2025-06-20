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
    WKBNamelist(; <keyword arguments>)

Construct a WKBNamelist instance, which holds parameters for the WKB ray tracing algorithm.

# Arguments:

  - `xrmin_dim::AbstractFloat = 0.0E+0`: Minimum x-coordinate of ray tracing domain in dimensional units.
  - `xrmax_dim::AbstractFloat = 1.0E+3`: Maximum x-coordinate of ray tracing domain in dimensional units.
  - `yrmin_dim::AbstractFloat = 0.0E+0`: Minimum y-coordinate of ray tracing domain in dimensional units.
  - `yrmax_dim::AbstractFloat = 1.0E+3`: Maximum y-coordinate of ray tracing domain in dimensional units.
  - `zrmin_dim::AbstractFloat = 0.0E+0`: Minimum z-coordinate of ray tracing domain in dimensional units.
  - `zrmax_dim::AbstractFloat = 1.0E+3`: Maximum z-coordinate of ray tracing domain in dimensional units.
  - `nrxl::Integer = 1`: Number of ray launch points in x-direction.
  - `nryl::Integer = 1`: Number of ray launch points in y-direction.
  - `nrzl::Integer = 1`: Number of ray launch points in z-direction.
  - `nrk_init::Integer = 1`: Initial number of wavenumber points in x-direction.
  - `nrl_init::Integer = 1`: Initial number of wavenumber points in y-direction.
  - `nrm_init::Integer = 1`: Initial number of wavenumber points in z-direction.
  - `nray_fac::Integer = 4`: Ray multiplication factor.
  - `fac_dk_init::AbstractFloat = 1.0E-1`: Initial wavenumber spacing factor in x-direction.
  - `fac_dl_init::AbstractFloat = 1.0E-1`: Initial wavenumber spacing factor in y-direction.
  - `fac_dm_init::AbstractFloat = 1.0E-1`: Initial wavenumber spacing factor in z-direction.
  - `branchr::Integer = -1`: Branch selection for ray tracing.
  - `merge_mode::AbstractMergeMode = ConstantWaveAction()`: Algorithm to use when merging ray beams.
  - `nsmth_wkb::Integer = 2`: Number of smoothing passes for WKB solution.
  - `lsmth_wkb::Bool = true`: Whether to apply smoothing to WKB solution.
  - `sm_filter::AbstractWKBFilter = Shapiro()`: Type of smoothing filter to apply to WKB solution.
  - `zmin_wkb_dim::AbstractFloat = 0.0`: Minimum height for WKB calculation in dimensional units.
  - `lsaturation::Bool = true`: Whether to apply wave saturation.
  - `alpha_sat::AbstractFloat = 1.0E+0`: Saturation coefficient.
  - `wkb_mode::AbstractWKBMode = MultiColumn()`: WKB calculation mode.
  - `blocking::Bool = false`: Whether to apply wave blocking.
  - `long_threshold::AbstractFloat = 2.5E-1`: Threshold for long wavelength approximation.
  - `drag_coefficient::AbstractFloat = 1.0E+0`: Coefficient for wave-induced drag.
  - `nwm::Integer = 1`: Number of wave modes.
  - `launch_algorithm::AbstractLaunchAlgorithm = Clip()`: Algorithm used for wave launching.
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
