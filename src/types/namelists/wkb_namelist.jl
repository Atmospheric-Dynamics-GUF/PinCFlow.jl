struct WKBNamelist{
  A <: AbstractFloat,
  B <: Integer,
  C <: AbstractMergeMode,
  D <: Bool,
  E <: AbstractWKBMode,
  F <: AbstractLaunchAlgorithm,
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
  cons_merge::C
  nsmth_wkb::B
  lsmth_wkb::D
  sm_filter::B
  zmin_wkb_dim::A
  lsaturation::D
  alpha_sat::A
  wkb_mode::E
  case_wkb::B
  blocking::D
  nwm::B
  launch_algorithm::F
end

function WKBNamelist(;
  xrmin_dim = 0.0,
  xrmax_dim = 1000.0,
  yrmin_dim = 0.0,
  yrmax_dim = 1000.0,
  zrmin_dim = 0.0,
  zrmax_dim = 1000.0,
  nrxl = 1,
  nryl = 1,
  nrzl = 1,
  nrk_init = 1,
  nrl_init = 1,
  nrm_init = 1,
  nray_fac = 4,
  fac_dk_init = 0.1,
  fac_dl_init = 0.1,
  fac_dm_init = 0.1,
  branchr = -1,
  cons_merge = ConstantWaveAction(),
  nsmth_wkb = 2,
  lsmth_wkb = true,
  sm_filter = 2,
  zmin_wkb_dim = 0.0,
  lsaturation = true,
  alpha_sat = 1.0,
  wkb_mode = MultiColumn(),
  case_wkb = 3,
  blocking = false,
  nwm = 1,
  launch_algorithm = Clip(),
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
    cons_merge,
    nsmth_wkb,
    lsmth_wkb,
    sm_filter,
    zmin_wkb_dim,
    lsaturation,
    alpha_sat,
    wkb_mode,
    case_wkb,
    blocking,
    nwm,
    launch_algorithm,
  )
end
