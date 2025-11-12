"""
```julia
WavePacketNamelist{A <: AbstractFloat, B <: Integer}
```

Namelist for the initialization of wave-packet test cases.

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
)::WavePacketNamelist
```

Construct a `WavePacketNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `wavepacketdim::B`: Dimensionality of the envelope (set to 1 for pure ``z``-dependence, 2 for ``x``- and ``z``-dependence, and 3 for dependence in all directions).

  - `lambdax_dim::A`: Zonal wavelength.

  - `lambday_dim::A`: Meridional wavelength.

  - `lambdaz_dim::A`: Vertical wavelength.

  - `x0_dim::A`: Zonal center.

  - `y0_dim::A`: Meridional center.

  - `z0_dim::A`: Vertical center.

  - `sigmax_dim::A`: Zonal envelope half-width.

  - `sigmay_dim::A`: Meridional envelope half-width.

  - `sigmaz_dim::A`: Vertical envelope half-width.

  - `a0::A`: Maximum amplitude in relation to saturation threshold.

  - `branch::B`: Frequency branch.
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
)::WavePacketNamelist
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
