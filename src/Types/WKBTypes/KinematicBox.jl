struct KinematicBox
    la::Vector{Int} # number of points in a (depends on kh)
    lq::Vector{Int} # number of points in q 
    lia::Vector{Float64} # logarithmic increment of the meshes in a
    liq::Vector{Float64} # logarithmic increment of the mesh in q
    loglia::Vector{Float64} # log(lia)
    logliq::Vector{Float64} # log(liq)
    aa::Vector{Vector{Float64}} # meshes of a
    qq::Vector{Vector{Float64}} # mesh of q
end

#dafualt initialization for NoWKB mode

function KinematicBox(wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad)::KinematicBox

    return KinematicBox(Int.(zeros(0)), Int.(zeros(0)), zeros(0), zeros(0), zeros(0), zeros(0), Vector{Vector{Float64}}(undef, 0), Vector{Vector{Float64}}(undef, 0))

end


function KinematicBox(amin::Vector{Float64}, 
    amax::Vector{Float64}, la::Vector{Int}, 
    qmin::Vector{Float64}, 
    qmax::Vector{Float64}, 
    lq::Vector{Int}, 
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{Triad2D, Triad3DIso})::KinematicBox
    @assert length(amin) == length(amax) && length(amin) == length(la) && length(qmin) == length(qmax) && length(qmin) == length(lq)
    mh=length(la)
    mq = length(lq)
    lia=Vector{Float64}(undef,mh)
    aa=Vector{Vector{Float64}}(undef,mh)
    liq=Vector{Float64}(undef,mq)
    qq=Vector{Vector{Float64}}(undef,mq)

    for ih in eachindex(la)
        aa[ih]=log_range(amin[ih], amax[ih], la[ih])
        lia[ih] = aa[ih][2]/aa[ih][1]
    end
    
    for iz in eachindex(lq)
        @assert lq[iz] >= 2 
            qq[iz] = log_range(qmin[iz], qmax[iz], lq[iz])
            liq[iz] = qq[iz][2] / qq[iz][1]
    end

    return KinematicBox(la, lq, lia, liq, log.(lia), log.(liq), aa, qq)

end
