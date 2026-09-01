function get_chi_cell_action end

function get_chi_cell_action(
    was_pred,
    delkp,
    delm,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    kpi::Integer,
    mi::Integer,
    x_size::Integer,
)
    was = was_pred[ii, jj, kk, kpi, mi]

    if !isfinite(was) || was <= 0.0
        return 0.0
    end

    spectral_measure =
        x_size == 1 ?
        abs(delm[mi]) :
        abs(delkp[kpi] * delm[mi])

    return was * spectral_measure
end