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
    xrmin_dim::AbstractFloat = 0.0E+0,
    xrmax_dim::AbstractFloat = 1.0E+3,
    yrmin_dim::AbstractFloat = 0.0E+0,
    yrmax_dim::AbstractFloat = 1.0E+3,
    zrmin_dim::AbstractFloat = 0.0E+0,
    zrmax_dim::AbstractFloat = 1.0E+3,
    nrxl::Integer = 1,
    nryl::Integer = 1,
    nrzl::Integer = 1,
    nrk_init::Integer = 1,
    nrl_init::Integer = 1,
    nrm_init::Integer = 1,
    nray_fac::Integer = 4,
    fac_dk_init::AbstractFloat = 1.0E-1,
    fac_dl_init::AbstractFloat = 1.0E-1,
    fac_dm_init::AbstractFloat = 1.0E-1,
    branchr::Integer = -1,
    merge_mode::AbstractMergeMode = ConstantWaveAction(),
    nsmth_wkb::Integer = 2,
    lsmth_wkb::Bool = true,
    sm_filter::AbstractWKBFilter = Shapiro(),
    zmin_wkb_dim::AbstractFloat = 0.0E+0,
    lsaturation::Bool = true,
    alpha_sat::AbstractFloat = 1.0E+0,
    wkb_mode::AbstractWKBMode = MultiColumn(),
    blocking::Bool = false,
    long_threshold::AbstractFloat = 2.5E-1,
    drag_coefficient::AbstractFloat = 1.0E+0,
    nwm::Integer = 1,
)
```

Construct a `WKBNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `xrmin_dim::A`: Lower bound in ``x`` of the region where ray volumes are initialized.

  - `xrmax_dim::A`: Upper bound in ``x`` of the region where ray volumes are initialized.

  - `yrmin_dim::A`: Lower bound in ``y`` of the region where ray volumes are initialized.

  - `yrmax_dim::A`: Upper bound in ``y`` of the region where ray volumes are initialized.

  - `zrmin_dim::A`: Lower bound in ``z`` of the region where ray volumes are initialized.

  - `zrmax_dim::A`: Upper bound in ``z`` of the region where ray volumes are initialized.

  - `nrxl::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{x}``-direction.

  - `nryl::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{y}``-direction.

  - `nrzl::B`: Number of ray-volumes launched per grid cell and wave mode in ``\\widehat{z}``-direction.

  - `nrk_init::B`: Number of ray-volumes launched per grid cell and wave mode in ``k``-direction.

  - `nrl_init::B`: Number of ray-volumes launched per grid cell and wave mode in ``l``-direction.

  - `nrm_init::B`: Number of ray-volumes launched per grid cell and wave mode in ``m``-direction.

  - `nray_fac::B`: Factor by which ray volumes are allowed to multiply in each dimension of physical space.

  - `fac_dk_init::A`: Relative inital ray-volume extent in ``k``.

  - `fac_dl_init::A`: Relative inital ray-volume extent in ``l``.

  - `fac_dm_init::A`: Relative inital ray-volume extent in ``m``.

  - `branchr::B`: Frequency branch.

  - `merge_mode::C`: Ray-volume merging strategy (conserved quantity).

  - `nsmth_wkb::B`: Order of the smoothing applied to the mean-flow tendencies.

  - `lsmth_wkb::D`: Switch for smoothing the mean-flow tendencies.

  - `sm_filter::E`: Filter to use for the smoothing of the mean-flow tendencies.

  - `zmin_wkb_dim::A`: Minimum altitude for ray-tracing and mean-flow impact.

  - `lsaturation::D`: Switch for the saturation scheme.

  - `alpha_sat::A`: Relative saturation threshold.

  - `wkb_mode::F`: Approximations used by MSGWaM.

  - `blocking::D`: Switch for parameterizing blocking in WKB-mountain-wave simulations.

  - `long_threshold::A`: Long-number threshold used by the blocked-layer scheme.

  - `drag_coefficient::A`: Dimensionless (relative) drag coefficient used by the blocked layer scheme.

  - `nwm::B`: Number of wave modes per grid cell.
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

function WKBNamelist(;
    xrmin_dim::AbstractFloat = 0.0E+0,
    xrmax_dim::AbstractFloat = 1.0E+3,
    yrmin_dim::AbstractFloat = 0.0E+0,
    yrmax_dim::AbstractFloat = 1.0E+3,
    zrmin_dim::AbstractFloat = 0.0E+0,
    zrmax_dim::AbstractFloat = 1.0E+3,
    nrxl::Integer = 1,
    nryl::Integer = 1,
    nrzl::Integer = 1,
    nrk_init::Integer = 1,
    nrl_init::Integer = 1,
    nrm_init::Integer = 1,
    nray_fac::Integer = 4,
    fac_dk_init::AbstractFloat = 1.0E-1,
    fac_dl_init::AbstractFloat = 1.0E-1,
    fac_dm_init::AbstractFloat = 1.0E-1,
    branchr::Integer = -1,
    merge_mode::AbstractMergeMode = ConstantWaveAction(),
    nsmth_wkb::Integer = 2,
    lsmth_wkb::Bool = true,
    sm_filter::AbstractWKBFilter = Shapiro(),
    zmin_wkb_dim::AbstractFloat = 0.0E+0,
    lsaturation::Bool = true,
    alpha_sat::AbstractFloat = 1.0E+0,
    wkb_mode::AbstractWKBMode = MultiColumn(),
    blocking::Bool = false,
    long_threshold::AbstractFloat = 2.5E-1,
    drag_coefficient::AbstractFloat = 1.0E+0,
    nwm::Integer = 1,
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
