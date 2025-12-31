function trapazoidal_with_logbin end

function trapazoidal_with_logbin(f::Vector{Float64},
     kk::Vector{Float64},
     lkk::Int, 
     lambda::Float64, 
     llambda::Float64, 
     imin::Int=2, 
     imax::Int=-1)::AbstractFloat
    f0 = 0.0
    inte = 0.0

    if imax == -1  #by default it integrates over the whole array kk
        imax = lkk
    end

    for i = imin:min(lkk, imax)
        if i == 1  # if imin = 1, to include the region [0, kk_min] using linear extrapolation
            f0 = max(0.0, -(f[2] - lambda * f[1]) / (lambda - 1.0)) #Linear extrapolation to k=0
            inte += 0.5 * kk[1] * (f[1] + f0)
        else
            inte += 0.5 * (f[i] * kk[i] + f[i-1] * kk[i-1]) * llambda
        end
    end

    return inte

end