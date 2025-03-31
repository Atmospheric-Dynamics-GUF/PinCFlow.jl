struct GridNamelist{A <: AbstractFloat, B <: Integer}
  mountainheight_dim::A
  mountainwidth_dim::A
  mountain_case::B
  range_factor::A
  spectral_modes::B
  envelope_reduction::A
  stretch_exponent::A
end

function GridNamelist(;
  mountainheight_dim = 100.0,
  mountainwidth_dim = 1000.0,
  mountain_case = 1,
  range_factor = 1.0,
  spectral_modes = 1,
  envelope_reduction = 0.0,
  stretch_exponent = 1.0,
)
  return GridNamelist(
    mountainheight_dim,
    mountainwidth_dim,
    mountain_case,
    range_factor,
    spectral_modes,
    envelope_reduction,
    stretch_exponent,
  )
end
