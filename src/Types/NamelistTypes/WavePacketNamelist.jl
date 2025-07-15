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
    u0_jet_dim::A 
    sigmaz_jet_dim::A
    z0_jet_dim::A
end

function WavePacketNamelist(;
    wavepacketdim = 1,
    lambdax_dim = 1.0E+3,
    lambday_dim = 0.0E+0,
    lambdaz_dim = 1.0E+3,
    x0_dim = 5.0E+3,
    y0_dim = 1.0E+3,
    z0_dim = 1.0E+4,
    sigmax_dim = 0.0E+0,
    sigmay_dim = 0.0E+0,
    sigmaz_dim = 2.0E+3,
    a0 = 9.0E-1,
    branch = -1,
    u0_jet_dim = 0.0E+0,
    sigmaz_jet_dim = 0.0E+0,
    z0_jet_dim = 0.0E+0,
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
        u0_jet_dim,
        sigmaz_jet_dim,
        z0_jet_dim,
    )
end
