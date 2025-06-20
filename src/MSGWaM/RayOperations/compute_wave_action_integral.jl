"""
    compute_wave_action_integral(merge_mode::ConstantWaveAction, nr::AbstractFloat, omegar::AbstractFloat, fxk::AbstractFloat, fyl::AbstractFloat, fzm::AbstractFloat)

# Arguments

  - `merge_mode::ConstantWaveAction`: Merge mode for constant wave action
  - `nr::AbstractFloat`: Wave action density
  - `omegar::AbstractFloat`: Intrinsic frequency
  - `fxk::AbstractFloat`: Phase space factor in x-k direction
  - `fyl::AbstractFloat`: Phase space factor in y-l direction
  - `fzm::AbstractFloat`: Phase space factor in z-m direction
"""
function compute_wave_action_integral(
    merge_mode::ConstantWaveAction,
    nr::AbstractFloat,
    omegar::AbstractFloat,
    fxk::AbstractFloat,
    fyl::AbstractFloat,
    fzm::AbstractFloat,
)
    return fxk * fyl * fzm * nr
end

"""
    compute_wave_action_integral(merge_mode::ConstantWaveEnergy, nr::AbstractFloat, omegar::AbstractFloat, fxk::AbstractFloat, fyl::AbstractFloat, fzm::AbstractFloat)

# Arguments

  - `merge_mode::ConstantWaveEnergy`: Merge mode for constant wave energy
  - `nr::AbstractFloat`: Wave action density
  - `omegar::AbstractFloat`: Intrinsic frequency
  - `fxk::AbstractFloat`: Phase space factor in x-k direction
  - `fyl::AbstractFloat`: Phase space factor in y-l direction
  - `fzm::AbstractFloat`: Phase space factor in z-m direction
"""
function compute_wave_action_integral(
    merge_mode::ConstantWaveEnergy,
    nr::AbstractFloat,
    omegar::AbstractFloat,
    fxk::AbstractFloat,
    fyl::AbstractFloat,
    fzm::AbstractFloat,
)
    return fxk * fyl * fzm * nr * omegar
end
