struct InterpCoef
    c_o::Array{Float64,2} # constant
    alphakp::Array{Float64,2} # -slope along horizontal 
    alpham::Array{Float64,2} # -slope along vertical 
    beta::Array{Float64,2} # nonlinear term coefficient
end

#dafualt initialization for NoWKB mode

function InterpCoef(wkb_mod::NoWKB)::InterpCoef
    
    v = Array{Float64}(undef, 0, 0)
    
    return InterpCoef(v, v, v, v)
end


function InterpCoef(kp::AbstractVector{<:AbstractFloat},
    m::AbstractVector{<:AbstractFloat}, 
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn})::InterpCoef

     v = zeros(length(kp), length(m))
    
    return InterpCoef(v, v, v, v)
end
