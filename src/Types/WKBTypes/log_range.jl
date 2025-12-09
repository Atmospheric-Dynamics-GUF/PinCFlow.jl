function log_range end

function log_range(kmin,kmax,M)
    logkmin = log(kmin)
    logkmax = log(kmax)
    logk = LinRange(logkmin,logkmax,M)
    k = exp.(logk)
    return k
end
