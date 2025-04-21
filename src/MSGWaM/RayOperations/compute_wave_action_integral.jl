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
