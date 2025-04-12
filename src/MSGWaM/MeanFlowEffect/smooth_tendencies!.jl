function smooth_tendencies!(state)
    (; sizex, sizey) = state.namelists.domain
    (; lsmth_wkb, sm_filter) = state.namelists.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.integrals

    if !lsmth_wkb
        return
    end

    if sizey == 1
        if sizex > 1
            smooth_tendencies!(dudt, state, sm_filter, XZ())
            smooth_tendencies!(dvdt, state, sm_filter, XZ())
            smooth_tendencies!(dthetadt, state, sm_filter, XZ())
        else
            error("Smoothing just in z not yet implemented.")
        end
    elseif sizex == 1
        smooth_tendencies!(dudt, state, sm_filter, YZ())
        smooth_tendencies!(dvdt, state, sm_filter, YZ())
        smooth_tendencies!(dthetadt, state, sm_filter, YZ())
    elseif sizex > 1
        smooth_tendencies!(dudt, state, sm_filter, XYZ())
        smooth_tendencies!(dvdt, state, sm_filter, XYZ())
        smooth_tendencies!(dthetadt, state, sm_filter, XYZ())
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Box, homogdir::XYZ)
    (; nx, ny, nz, nbx, nby, nbz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    if min(nx, nbx) < nsmth
        error("min(nx, nbx) too small for smoothing")
    end
    if min(ny, nby) < nsmth
        error("min(ny, nby) too small for smoothing")
    end
    if min(nz, nbz) < nsmth
        error("min(nz, nbz) too small for smoothing")
    end

    flxwkb_0 = copy(flxwkb)

    for k in k0:k1, j in j0:j1, i in i0:i1
        flxwkb[i, j, k] =
            sum(
                flxwkb_0[
                    (i - nsmth):(i + nsmth),
                    (j - nsmth):(j + nsmth),
                    (k - nsmth):(k + nsmth),
                ],
            ) / real((2 * nsmth + 1)^3)
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Box, homogdir::XZ)
    (; nx, ny, nz, nbx, nby, nbz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    if min(nx, nbx) < nsmth
        error("min(nx, nbx) too small for smoothing")
    end
    if min(nz, nbz) < nsmth
        error("min(nz, nbz) too small for smoothing")
    end

    flxwkb_0 = copy(flxwkb)

    for k in k0:k1, j in j0:j1, i in i0:i1
        flxwkb[i, j, k] =
            sum(flxwkb_0[(i - nsmth):(i + nsmth), j, (k - nsmth):(k + nsmth)]) /
            real((2 * nsmth + 1)^2)
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Box, homogdir::YZ)
    (; nx, ny, nz, nbx, nby, nbz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    if min(ny, nby) < nsmth
        error("min(ny, nby) too small for smoothing")
    end
    if min(nz, nbz) < nsmth
        error("min(nz, nbz) too small for smoothing")
    end

    flxwkb_0 = copy(flxwkb)

    for k in k0:k1, j in j0:j1, i in i0:i1
        flxwkb[i, j, k] =
            sum(flxwkb_0[i, (j - nsmth):(j + nsmth), (k - nsmth):(k + nsmth)]) /
            real((2 * nsmth + 1)^2)
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Box, homogdir::X)
    (; nx, ny, nz, nbx, nby, nbz, i0, i1, j0, j1, k0, k1) = state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    flxwkb_0 = copy(flxwkb)

    for k in k0:k1, j in j0:j1, i in i0:i1
        flxwkb[i, j, k] = sum(flxwkb_0[i0:i1, j, k]) / real(nx)
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Shapiro, homogdir::XYZ)
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    if min(nx, nbx) < nsmth
        error("min(nx, nbx) too small for smoothing")
    end
    if min(ny, nby) < nsmth
        error("min(ny, nby) too small for smoothing")
    end
    if min(nz, nbz) < nsmth
        error("min(nz, nbz) too small for smoothing")
    end

    flxwkb_0 = copy(flxwkb)
    flxwkb_1 = zeros(flxwkb)

    if nsmth == 1
        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in (j0 - nby):(j1 + nby), i in i0:i1
            flxwkb_0[i, j, k] =
                (
                    flxwkb[i - 1, j, k] +
                    flxwkb[i + 1, j, k] +
                    2.0 * flxwkb[i, j, k]
                ) / 4.0
        end

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    flxwkb_0[i, j - 1, k] +
                    flxwkb_0[i, j + 1, k] +
                    2.0 * flxwkb_0[i, j, k]
                ) / 4.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 1] +
                    flxwkb_1[i, j, k + 1] +
                    2.0 * flxwkb_1[i, j, k]
                ) / 4.0
        end

    elseif nsmth == 2

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in (j0 - nby):(j1 + nby), i in i0:i1
            flxwkb_0[i, j, k] =
                (
                    -flxwkb[i - 2, j, k] - flxwkb[i + 2, j, k] +
                    4.0 * (flxwkb[i - 1, j, k] + flxwkb[i + 1, j, k]) +
                    10.0 * flxwkb[i, j, k]
                ) / 16.0
        end

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    -flxwkb_0[i, j - 2, k] - flxwkb_0[i, j + 2, k] +
                    4.0 * (flxwkb_0[i, j - 1, k] + flxwkb_0[i, j + 1, k]) +
                    10.0 * flxwkb_0[i, j, k]
                ) / 16.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    -flxwkb_1[i, j, k - 2] - flxwkb_1[i, j, k + 2] +
                    4.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    10.0 * flxwkb_1[i, j, k]
                ) / 16.0
        end

    elseif nsmth == 3

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in (j0 - nby):(j1 + nby), i in i0:i1
            flxwkb_0[i, j, k] =
                (
                    flxwkb[i - 3, j, k] + flxwkb[i + 3, j, k] -
                    6.0 * (flxwkb[i - 2, j, k] + flxwkb[i + 2, j, k]) +
                    15.0 * (flxwkb[i - 1, j, k] + flxwkb[i + 1, j, k]) +
                    44.0 * flxwkb[i, j, k]
                ) / 64.0
        end

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    flxwkb_0[i, j - 3, k] + flxwkb_0[i, j + 3, k] -
                    6.0 * (flxwkb_0[i, j - 2, k] + flxwkb_0[i, j + 2, k]) +
                    15.0 * (flxwkb_0[i, j - 1, k] + flxwkb_0[i, j + 1, k]) +
                    44.0 * flxwkb_0[i, j, k]
                ) / 64.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 3] + flxwkb_1[i, j, k + 3] -
                    6.0 * (flxwkb_1[i, j, k - 2] + flxwkb_1[i, j, k + 2]) +
                    15.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    44.0 * flxwkb_1[i, j, k]
                ) / 64.0
        end

    elseif nsmth == 4

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in (j0 - nby):(j1 + nby), i in i0:i1
            flxwkb_0[i, j, k] =
                (
                    -flxwkb[i - 4, j, k] - flxwkb[i + 4, j, k] +
                    8.0 * (flxwkb[i - 3, j, k] + flxwkb[i + 3, j, k]) -
                    28.0 * (flxwkb[i - 2, j, k] + flxwkb[i + 2, j, k]) +
                    56.0 * (flxwkb[i - 1, j, k] + flxwkb[i + 1, j, k]) +
                    186.0 * flxwkb[i, j, k]
                ) / 256.0
        end

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    -flxwkb_0[i, j - 4, k] - flxwkb_0[i, j + 4, k] +
                    8.0 * (flxwkb_0[i, j - 3, k] + flxwkb_0[i, j + 3, k]) -
                    28.0 * (flxwkb_0[i, j - 2, k] + flxwkb_0[i, j + 2, k]) +
                    56.0 * (flxwkb_0[i, j - 1, k] + flxwkb_0[i, j + 1, k]) +
                    186.0 * flxwkb_0[i, j, k]
                ) / 256.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 4] +
                    flxwkb_1[i, j, k + 4] +
                    8.0 * (flxwkb_1[i, j, k - 3] + flxwkb_1[i, j, k + 3]) -
                    28.0 * (flxwkb_1[i, j, k - 2] + flxwkb_1[i, j, k + 2]) +
                    56.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    186.0 * flxwkb_1[i, j, k]
                ) / 256.0
        end

    else
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                sum(
                    flxwkb_0[
                        (i - nsmth):(i + nsmth),
                        (j - nsmth):(j + nsmth),
                        (k - nsmth):(k + nsmth),
                    ],
                ) / real((2 * nsmth + 1)^3)
        end
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Shapiro, homogdir::XZ)
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    if min(nx, nbx) < nsmth
        error("min(nx, nbx) too small for smoothing")
    end
    if min(nz, nbz) < nsmth
        error("min(nz, nbz) too small for smoothing")
    end

    flxwkb_0 = copy(flxwkb)
    flxwkb_1 = zeros(flxwkb)

    if nsmth == 1

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    flxwkb_0[i - 1, j, k] +
                    flxwkb_0[i + 1, j, k] +
                    2.0 * flxwkb_0[i, j, k]
                ) / 4.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 1] +
                    flxwkb_1[i, j, k + 1] +
                    2.0 * flxwkb_1[i, j, k]
                ) / 4.0
        end

    elseif nsmth == 2

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    -flxwkb_0[i - 2, j, k] - flxwkb_0[i + 2, j, k] +
                    4.0 * (flxwkb_0[i - 1, j, k] + flxwkb_0[i + 1, j, k]) +
                    10.0 * flxwkb_0[i, j, k]
                ) / 16.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    -flxwkb_1[i, j, k - 2] - flxwkb_1[i, j, k + 2] +
                    4.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    10.0 * flxwkb_1[i, j, k]
                ) / 16.0
        end

    elseif nsmth == 3

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    flxwkb_0[i - 3, j, k] + flxwkb_0[i + 3, j, k] -
                    6.0 * (flxwkb_0[i - 2, j, k] + flxwkb_0[i + 2, j, k]) +
                    15.0 * (flxwkb_0[i - 1, j, k] + flxwkb_0[i + 1, j, k]) +
                    44.0 * flxwkb_0[i, j, k]
                ) / 64.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 3] + flxwkb_1[i, j, k + 3] -
                    6.0 * (flxwkb_1[i, j, k - 2] + flxwkb_1[i, j, k + 2]) +
                    15.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    44.0 * flxwkb_1[i, j, k]
                ) / 64.0
        end

    elseif nsmth == 4

        # smooth in x
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    -flxwkb_0[i - 4, j, k] - flxwkb_0[i + 4, j, k] +
                    8.0 * (flxwkb_0[i - 3, j, k] + flxwkb_0[i + 3, j, k]) -
                    28.0 * (flxwkb_0[i - 2, j, k] + flxwkb_0[i + 2, j, k]) +
                    56.0 * (flxwkb_0[i - 1, j, k] + flxwkb_0[i + 1, j, k]) +
                    186.0 * flxwkb_0[i, j, k]
                ) / 256.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 4] +
                    flxwkb_1[i, j, k + 4] +
                    8.0 * (flxwkb_1[i, j, k - 3] + flxwkb_1[i, j, k + 3]) -
                    28.0 * (flxwkb_1[i, j, k - 2] + flxwkb_1[i, j, k + 2]) +
                    56.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    186.0 * flxwkb_1[i, j, k]
                ) / 256.0
        end

    else
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                sum(
                    flxwkb_0[
                        (i - nsmth):(i + nsmth),
                        j,
                        (k - nsmth):(k + nsmth),
                    ],
                ) / real((2 * nsmth + 1)^2)
        end
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Shapiro, homogdir::YZ)
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    if min(ny, nby) < nsmth
        error("min(ny, nby) too small for smoothing")
    end
    if min(nz, nbz) < nsmth
        error("min(nz, nbz) too small for smoothing")
    end

    flxwkb_0 = copy(flxwkb)
    flxwkb_1 = zeros(flxwkb)

    if nsmth == 1

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    flxwkb_0[i, j - 1, k] +
                    flxwkb_0[i, j + 1, k] +
                    2.0 * flxwkb_0[i, j, k]
                ) / 4.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 1] +
                    flxwkb_1[i, j, k + 1] +
                    2.0 * flxwkb_1[i, j, k]
                ) / 4.0
        end

    elseif nsmth == 2

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    -flxwkb_0[i, j - 2, k] - flxwkb_0[i, j + 2, k] +
                    4.0 * (flxwkb_0[i, j - 1, k] + flxwkb_0[i, j + 1, k]) +
                    10.0 * flxwkb_0[i, j, k]
                ) / 16.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    -flxwkb_1[i, j, k - 2] - flxwkb_1[i, j, k + 2] +
                    4.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    10.0 * flxwkb_1[i, j, k]
                ) / 16.0
        end

    elseif nsmth == 3

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    flxwkb_0[i, j - 3, k] + flxwkb_0[i, j + 3, k] -
                    6.0 * (flxwkb_0[i, j - 2, k] + flxwkb_0[i, j + 2, k]) +
                    15.0 * (flxwkb_0[i, j - 1, k] + flxwkb_0[i, j + 1, k]) +
                    44.0 * flxwkb_0[i, j, k]
                ) / 64.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 3] + flxwkb_1[i, j, k + 3] -
                    6.0 * (flxwkb_1[i, j, k - 2] + flxwkb_1[i, j, k + 2]) +
                    15.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    44.0 * flxwkb_1[i, j, k]
                ) / 64.0
        end

    elseif nsmth == 4

        # smooth in y
        for k in (k0 - nbz):(k1 + nbz), j in j0:j1, i in i0:i1
            flxwkb_1[i, j, k] =
                (
                    -flxwkb_0[i, j - 4, k] - flxwkb_0[i, j + 4, k] +
                    8.0 * (flxwkb_0[i, j - 3, k] + flxwkb_0[i, j + 3, k]) -
                    28.0 * (flxwkb_0[i, j - 2, k] + flxwkb_0[i, j + 2, k]) +
                    56.0 * (flxwkb_0[i, j - 1, k] + flxwkb_0[i, j + 1, k]) +
                    186.0 * flxwkb_0[i, j, k]
                ) / 256.0
        end

        # smooth in z
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                (
                    flxwkb_1[i, j, k - 4] +
                    flxwkb_1[i, j, k + 4] +
                    8.0 * (flxwkb_1[i, j, k - 3] + flxwkb_1[i, j, k + 3]) -
                    28.0 * (flxwkb_1[i, j, k - 2] + flxwkb_1[i, j, k + 2]) +
                    56.0 * (flxwkb_1[i, j, k - 1] + flxwkb_1[i, j, k + 1]) +
                    186.0 * flxwkb_1[i, j, k]
                ) / 256.0
        end

    else
        for k in k0:k1, j in j0:j1, i in i0:i1
            flxwkb[i, j, k] =
                sum(
                    flxwkb_0[
                        i,
                        (j - nsmth):(j + nsmth),
                        (k - nsmth):(k + nsmth),
                    ],
                ) / real((2 * nsmth + 1)^2)
        end
    end
    return
end

function smooth_tendencies!(flxwkb, state, sm_filter::Shapiro, homogdir::X)
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; nsmth_wkb) = state.namelists.WKBNamelist

    nsmth = nsmth_wkb

    flxwkb_0 = copy(flxwkb)

    for k in k0:k1, j in j0:j1, i in i0:i1
        flxwkb[i, j, k] = sum(flxwkb_0[i0:i1, j, k]) / real(nx)
    end
    return
end
