function get_nl_time_scale! end

function get_nl_time_scale!(spec_tend::TriadTendencies,
    ii::Integer,
    jj::Integer,
    kk::Integer,)::AbstractFloat
    inds = spec_tend.wavespectrum[ii, jj, kk, :, :] .> 1.e-40
    spec_tend.diag_time .= 0
    if any(inds)
        @. spec_tend.diag_time[inds] = spec_tend.col_int[inds] ./ spec_tend.wavespectrum[ii, jj, kk, :, :][inds]
        return 1.0 / maximum(abs.(spec_tend.diag_time[inds]))
    else
        return -1.0
    end
end