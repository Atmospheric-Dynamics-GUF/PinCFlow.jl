struct InterpCoef
    c_o::Array{Float64,2} # constant
    alphakp::Array{Float64,2} # -slope along horizontal 
    alpham::Array{Float64,2} # -slope along vertical 
    beta::Array{Float64,2} # nonlinear term coefficient
end

#dafualt initialization for NoWKB mode

function InterpCoef(wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
   triad_mode::NoTriad)::InterpCoef
    
     Array{Float64}(undef, 0, 0)
    
    return InterpCoef(Array{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0), Array{Float64}(undef, 0, 0))
end


function InterpCoef(kk::AbstractVector{<:AbstractFloat},
    mm::AbstractVector{<:AbstractFloat}, 
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{Triad2D, Triad3DIso})::InterpCoef

     c_o = zeros(length(kk), length(mm))
     alphakp = zeros(length(kk), length(mm))
     alpham = zeros(length(kk), length(mm))
     beta = zeros(length(kk), length(mm))
    
    return InterpCoef(c_o, alphakp, alpham, beta)
end
