"""
```julia
WavePacketNamelist{A <: AbstractFloat, B <: Integer}
```
"""
struct WavePacketNamelist{A <: AbstractFloat, B <: Integer}
    wavepacketdim::B
    lambdax_dim::A
    lambday_dim::A
    lambdaz_dim::A
    x0_dim::A
    y0_dim::A
    z0_dim::A
    sigmax_dim::A
    sigmay_dim::A
    sigmaz_dim::A
    a0::A
    branch::B
end

"""
```julia
WavePacketNamelist(;
    wavepacketdim::Integer = 1,
    lambdax_dim::AbstractFloat = 1.0E+3,
    lambday_dim::AbstractFloat = 0.0E+0,
    lambdaz_dim::AbstractFloat = 1.0E+3,
    x0_dim::AbstractFloat = 5.0E+3,
    y0_dim::AbstractFloat = 1.0E+3,
    z0_dim::AbstractFloat = 1.0E+4,
    sigmax_dim::AbstractFloat = 0.0E+0,
    sigmay_dim::AbstractFloat = 0.0E+0,
    sigmaz_dim::AbstractFloat = 2.0E+3,
    a0::AbstractFloat = 9.0E-1,
    branch::Integer = -1,
)
```
"""
function WavePacketNamelist(;
    wavepacketdim::Integer = 1,
    lambdax_dim::AbstractFloat = 1.0E+3,
    lambday_dim::AbstractFloat = 0.0E+0,
    lambdaz_dim::AbstractFloat = 1.0E+3,
    x0_dim::AbstractFloat = 5.0E+3,
    y0_dim::AbstractFloat = 1.0E+3,
    z0_dim::AbstractFloat = 1.0E+4,
    sigmax_dim::AbstractFloat = 0.0E+0,
    sigmay_dim::AbstractFloat = 0.0E+0,
    sigmaz_dim::AbstractFloat = 2.0E+3,
    a0::AbstractFloat = 9.0E-1,
    branch::Integer = -1,
)
    return WavePacketNamelist(
        wavepacketdim,
        lambdax_dim,
        lambday_dim,
        lambdaz_dim,
        x0_dim,
        y0_dim,
        z0_dim,
        sigmax_dim,
        sigmay_dim,
        sigmaz_dim,
        a0,
        branch,
    )
end
