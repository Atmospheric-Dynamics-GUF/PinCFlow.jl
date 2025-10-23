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
struct MultiWavePacketNamelist{A <: AbstractArray{<:AbstractFloat, 1}, B <: AbstractArray{<:Integer, 1}, C <: Integer } 
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
    nwm::C
end

function MultiWavePacketNamelist(;
    nwm = 1,
    wavepacketdim = [1],
    lambdax_dim = [1.0E+3],
    lambday_dim = [0.0E+0],
    lambdaz_dim = [1.0E+3],
    x0_dim = [5.0E+3],
    y0_dim = [1.0E+3],
    z0_dim = [1.0E+4],
    sigmax_dim = [0.0E+0],
    sigmay_dim = [0.0E+0],
    sigmaz_dim = [2.0E+3],
    a0 = [9.0E-1],
)::MultiWavePacketNamelist
    return MultiWavePacketNamelist(
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
        nwm,
    )
end

function MultiWavePacketNamelist(;random_wavepackets::RandomWavePackets, number_wavepackets::Integer, domain::DomainNamelist, testcase::AbstractTestCase)

    nwm = number_wavepackets
     wavepacketdim = ones(Int, nwm)

    # define range for random wavelengths x-direction
    int_b_lambdax_dim = 1.e4
    int_e_lambdax_dim = 1.e5

    if testcase isa MultipleWavePackets
        #check if resolution is sufficient
        if ((lx_dim[2]-lx_dim[1])/sizex > int_b_lambdax_dim/10.)
            println("Resolution in x-direction too coarse for random wavepackets in MultipleWavePackets testcase !!! ")
            exit(1)
        end
        if ((lz_dim[2]-lz_dim[1])/sizez > int_b_lambdaz_dim/10.)
            println("Resolution in z-direction too coarse for random wavepackets in MultipleWavePackets testcase !!! ")
            exit(1)
        end    

    elseif !(testcase isa WKBMultipleWavePackets)
        println("only testcase=MultipleWavePackets or WKBMultipleWavePacket for  random wavepackets !!! ")
        exit(1)    

    end

    # define range for random wavelengths z-direction
    int_b_lambdaz_dim = 1.e3
    int_e_lambdaz_dim = 2.e3

    int_b_x0_dim = domain.lx_dim[2]/2.
    int_e_x0_dim = domain.lx_dim[2]/2.

    z0_issr = 8.e3 # [m] to be consistent with the value in IcePredictands.jl
    int_b_z0_dim = z0_issr - 2.e+3
    int_e_z0_dim = z0_issr + 2.e+3
    
    # random width envelope
    int_b_sigmaz_dim = 1.0E+4
    int_e_sigmaz_dim = 2.0E+4

    # random amplitude
    int_b_a0 = 0.1E+0
    int_e_a0 = 1.0E+0

    a0 = int_b_a0 .+ (int_e_a0 - int_b_a0) .* rand(Float64, nwm)

    x0_dim = int_b_x0_dim .+ (int_e_x0_dim - int_b_x0_dim) .* rand(Float64, nwm)
    y0_dim = 1.e+3 .* ones(Float64, nwm)
    z0_dim = int_b_z0_dim .+ (int_e_z0_dim - int_b_z0_dim) .* rand(Float64, nwm)

    lambdax_dim = int_b_lambdax_dim .+ (int_e_lambdax_dim - int_b_lambdax_dim) .* rand(Float64, nwm)
    lambday_dim = zeros(Float64, nwm)
    lambdaz_dim = int_b_lambdaz_dim .+ (int_e_lambdaz_dim - int_b_lambdaz_dim) .* rand(Float64, nwm)

     sigmax_dim = zeros(Float64, nwm)
     sigmay_dim = zeros(Float64, nwm)
     sigmaz_dim = int_b_sigmaz_dim .+ (int_e_sigmaz_dim - int_b_sigmaz_dim) .* rand(Float64, nwm)

    return MultiWavePacketNamelist(
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
        nwm,
    )
end    

