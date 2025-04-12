function merge_wave_action(
    merge_mode::ConstantWaveAction,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    nr::AbstractFloat,
    omegar::AbstractFloat,
)
    return axk * ayl * azm * nr
end

function merge_wave_action(
    merge_mode::ConstantWaveEnergy,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    nr::AbstractFloat,
    omegar::AbstractFloat,
)
    return axk * ayl * azm * nr * omegar
end

function merge_wave_action(
    merge_mode::ConstantWaveAction,
    self::MergedRays,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    nr::AbstractFloat,
    omegar::AbstractFloat,
)
    return self.nr + axk * ayl * azm * nr
end

function merge_wave_action(
    merge_mode::ConstantWaveAction,
    self::MergedRays,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    nr::AbstractFloat,
    omegar::AbstractFloat,
)
    return self.nr + axk * ayl * azm * nr * omegar
end
