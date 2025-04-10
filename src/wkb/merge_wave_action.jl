function merge_wave_action(
    merge_mode::ConstantWaveAction,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    dens::AbstractFloat,
    omir::AbstractFloat,
)
    return axk * ayl * azm * dens
end

function merge_wave_action(
    merge_mode::ConstantWaveEnergy,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    dens::AbstractFloat,
    omir::AbstractFloat,
)
    return axk * ayl * azm * dens * omir
end

function merge_wave_action(
    merge_mode::ConstantWaveAction,
    self::MergedRayVolume,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    dens::AbstractFloat,
    omir::AbstractFloat,
)
    return self.dens + axk * ayl * azm * dens
end

function merge_wave_action(
    merge_mode::ConstantWaveAction,
    self::MergedRayVolume,
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    dens::AbstractFloat,
    omir::AbstractFloat,
)
    return self.dens + axk * ayl * azm * dens * omir
end
