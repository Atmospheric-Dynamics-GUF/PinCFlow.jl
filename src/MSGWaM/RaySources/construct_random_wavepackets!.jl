function construct_random_wavepackets!(wavepacket_namelist::MultiWavePacketNamelist, 
domain::DomainNamelist, testcase::AbstractTestCase)

    (; lx_dim, ly_dim, lz_dim) = domain
    (; sizex, sizey, sizez) = domain

    random_wavepackets = wavepacket_namelist.random_wavepackets
    nwm = wavepacket_namelist.nwm
    wavepacketdim = ones(Int, nwm)

    # define range for random wavelengths x-direction
    int_b_lambdax_dim = 1.e4
    int_e_lambdax_dim = 1.e5

    # define range for random wavelengths y-direction
    int_b_lambday_dim = 1.e4
    int_e_lambday_dim = 1.e5

    # define range for random wavelengths z-direction
    int_b_lambdaz_dim = 1.e3
    int_e_lambdaz_dim = 1.e4

    if testcase isa MultipleWavePackets

        #check if resolution is sufficient for resolving the waves
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

    # define intervals for centers wavepackets
    int_b_x0_dim = lx_dim[2]/2.
    int_e_x0_dim = lx_dim[2]/2.

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
    # random sign
    lambdaz_dim .*= 2.0 .*rand(Bool, nwm) .- 1

    sigmax_dim = zeros(Float64, nwm)
    sigmay_dim = zeros(Float64, nwm)
    sigmaz_dim = int_b_sigmaz_dim .+ (int_e_sigmaz_dim - int_b_sigmaz_dim) .* rand(Float64, nwm)

    (; wavepacketdim,
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
        random_wavepackets,
    ) = wavepacket_namelist 
  
    return 
end    

