struct ElasticModeSelection{
    A <: AbstractVector{<:Integer},
    B <: AbstractMatrix{<:Integer},
    C <: AbstractMatrix{<:AbstractFloat},
}
    sorted_wave_mode_indices::A
    launch_mode_count::B
    launch_power_fraction::C
end

function ElasticModeSelection(
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
)::ElasticModeSelection
    return ElasticModeSelection(
        zeros(Int, wave_modes),
        zeros(Int, nxx, nyy),
        zeros(nxx, nyy),
    )
end
