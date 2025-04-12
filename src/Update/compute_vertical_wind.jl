function compute_vertical_wind(
    i::Integer,
    j::Integer,
    k::Integer,
    predictands::Predictands,
    grid::Grid,
)
    (; u, v, w) = predictands

    uedger = u[i, j, k]
    uuedger = u[i, j, k + 1]
    uedgel = u[i - 1, j, k]
    uuedgel = u[i - 1, j, k + 1]
    vedgef = v[i, j, k]
    vuedgef = v[i, j, k + 1]
    vedgeb = v[i, j - 1, k]
    vuedgeb = v[i, j - 1, k + 1]
    wedgeu = w[i, j, k]

    return transform(
        i,
        j,
        k,
        uedger,
        uuedger,
        uedgel,
        uuedgel,
        vedgef,
        vuedgef,
        vedgeb,
        vuedgeb,
        wedgeu,
        Cartesian(),
        grid,
    )
end
