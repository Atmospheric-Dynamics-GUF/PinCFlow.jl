"""
```julia
compute_wave_action_integral(
    merge_mode::ConstantWaveAction,
    nr::AbstractFloat,
    omegar::AbstractFloat,
    fxk::AbstractFloat,
    fyl::AbstractFloat,
    fzm::AbstractFloat,
)
```

Returns the wave action obtained by multiplying the given phase-space wave-action density with the given phase-space volume.

This method is used to implement conservation of wave action in ray-volume merging.

```julia
compute_wave_action_integral(
    merge_mode::ConstantWaveEnergy,
    nr::AbstractFloat,
    omegar::AbstractFloat,
    fxk::AbstractFloat,
    fyl::AbstractFloat,
    fzm::AbstractFloat,
)
```

Returns the wave energy obtained by multiplying the given phase-space wave-action density with the given intrinsic frequency and phase-space volume.

This method is used to implement conservation of wave energy in ray-volume merging.

# Arguments

  - `merge_mode`: Merging strategy.

  - `nr`: Phase-space wave-action density.

  - `omegar`: Intrinsic frequency.

  - `fxk`: Phase space factor in the ``x``-``k`` subspace.

  - `fyl`: Phase space factor in the ``y``-``l`` subspace.

  - `fzm`: Phase space factor in the ``z``-``m`` subspace.

# Returns

  - `::AbstractFloat`: Either `fxk * fyl * fzm * nr` (wave action) or `fxk * fyl * fzm * nr * omegar` (wave energy), depending on the method.
"""
function compute_wave_action_integral end

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
