struct KinematicBox
    la::Vector{Int} # number of points in a (depends on kh)
    lq::Int # number of points in q 
    lia::Vector{Float64} # logarithmic increment of the meshes in a
    liq::Float64 # logarithmic increment of the mesh in q
    loglia::Vector{Float64} # log(lia)
    logliq::Float64 # log(liq)
    aa::Vector{Vector{Float64}} # meshes of a
    qq::Vector{Float64} # mesh of q
end

#dafualt initialization for NoWKB mode

function KinematicBox(wkb_mod::NoWKB)::KinematicBox
    
    v = Vector{Vector{Float64}}(undef, 0)
    
    return KinematicBox(Int.(zeros(0)), 0, zeros(0), 0.0, zeros(0), 0.0, v, zeros(0))
end


function KinematicBox(amin::Vector{Float64}, amax::Vector{Float64}, la::Vector{Int}, qmin::Float64, qmax::Float64, lq::Int, 
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn})::KinematicBox
    @assert length(amin) == length(amax) && length(amin) == length(la)
    mh=length(la)
    lia=Vector{Float64}(undef,mh)
    aa=Vector{Vector{Float64}}(undef,mh)

    for ih in eachindex(la)
        aa[ih]=log_range(amin[ih],amax[ih],la[ih])
        lia[ih] = aa[ih][2]/aa[ih][1]
    end
    
    qq=log_range(qmin, qmax, lq) 
    liq=qq[2]/qq[1]

    return KinematicBox(la, lq, lia, liq, log.(lia), log(liq), aa, qq)
end
