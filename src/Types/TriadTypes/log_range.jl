function log_range end

function log_range(kmin,kmax,M)
    @assert (kmin > 0) && (kmax >=0) && M >= 2
    if kmax == 0
        k = [0.0]
    elseif kmax <= kmin
        k = [kmin]
    else
        logkmin = log(kmin)
        logkmax = log(kmax)
        logk = LinRange(logkmin,logkmax,M)
        k = exp.(logk)
    end
    return k
end
