abstract type AbstractMergeMode end
struct ConstantWaveEnergy <: AbstractMergeMode end
struct ConstantWaveAction <: AbstractMergeMode end

abstract type AbstractLaunchAlgorithm end
struct Clip <: AbstractLaunchAlgorithm end
struct Scale <: AbstractLaunchAlgorithm end

struct WKBNamelist{
  A <: AbstractFloat,
  B <: Integer,
  C <: AbstractMergeMode,
  D <: Bool,
  E <: AbstractLaunchAlgorithm,
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
  single_column::D
  steady_state::D
  case_wkb::B
  amp_wkb::A
  wlrx_init::A
  wlry_init::A
  wlrz_init::A
  xr0_dim::A
  yr0_dim::A
  zr0_dim::A
  sigwpx_dim::A
  sigwpy_dim::A
  sigwpz_dim::A
  blocking::D
  nwm::B
  launch_algorithm::E
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
  single_column = false,
  steady_state = false,
  case_wkb = 3,
  amp_wkb = 1.0,
  wlrx_init = 100.0,
  wlry_init = 100.0,
  wlrz_init = 100.0,
  xr0_dim = 500.0,
  yr0_dim = 500.0,
  zr0_dim = 500.0,
  sigwpx_dim = 500.0,
  sigwpy_dim = 500.0,
  sigwpz_dim = 500.0,
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
    single_column,
    steady_state,
    case_wkb,
    amp_wkb,
    wlrx_init,
    wlry_init,
    wlrz_init,
    xr0_dim,
    yr0_dim,
    zr0_dim,
    sigwpx_dim,
    sigwpy_dim,
    sigwpz_dim,
    blocking,
    nwm,
    launch_algorithm,
  )
end
