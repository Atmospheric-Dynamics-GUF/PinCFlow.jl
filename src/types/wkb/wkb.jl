struct WKB{
    A <: Integer,
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: Rays,
    D <: SurfaceIndices,
    E <: Increments,
    F <: Integrals,
    G <: AbstractMatrix{<:AbstractFloat},
    H <: Forces,
    I <: AbstractFloat,
}
    nxray::A
    nyray::A
    nzray::A
    nxray_wrk::A
    nyray_wrk::A
    nzray_wrk::A
    nray_max::A
    nray_wrk::A
    n_sfc::A
    nray::B
    rays::C
    surface_indices::D
    increments::E
    integrals::F
    cgx_max::I
    cgy_max::I
    cgz_max::B
    zb::G
    gwmomforce::H
end

function WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    atmosphere::Atmosphere,
    predictands::Predictands,
)
    (; testcase) = namelists.setting
    return WKB(
        namelists,
        constants,
        domain,
        grid,
        atmosphere,
        predictands,
        testcase,
    )
end

function WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    atmosphere::Atmosphere,
    predictands::Predictands,
    testcase::AbstractTestCase,
)
    return WKB(
        [0 for i in 1:9]...,
        zeros(0, 0, 0),
        Rays(0, 0, 0, 0),
        SurfaceIndices(0, 0, 0),
        Increments(0, 0, 0, 0),
        Integrals(0, 0, 0),
        zeros(0, 0, 0),
        Forces(0, 0, 0),
    )
end

function WKB(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    atmosphere::Atmosphere,
    predictands::Predictands,
    testcase::AbstractWKBTestCase,
)

    # Get necessary fields.
    (;
        xrmin_dim,
        xrmax_dim,
        yrmin_dim,
        yrmax_dim,
        zrmin_dim,
        zrmax_dim,
        nrxl,
        nryl,
        nrzl,
        nrk_init,
        nrl_init,
        nrm_init,
        nray_fac,
        nwm,
        wkb_mode,
        fac_dk_init,
        fac_dl_init,
        fac_dm_init,
    ) = namelists.wkb
    (; tref, lref) = constants
    (; sizex, sizey, sizez) = namelists.domain
    (; f_coriolis_dim) = namelists.atmosphere
    (; nxx, nyy, nzz, io, jo, i0, i1, j0, j1, k0, k1, master, comm) = domain
    (; lx, ly, lz, dx, dy) = grid

    f_cor_nd = f_coriolis_dim * tref
    # Set Coriolis parameter.

    # Set boundaries for ray-volume propagation.
    xrmin = xrmin_dim / lref
    xrmax = xrmax_dim / lref
    yrmin = yrmin_dim / lref
    yrmax = yrmax_dim / lref
    zrmin = zrmin_dim / lref
    zrmax = zrmax_dim / lref

    # Check if the boundaries for ray-volume propagation are within the domain.
    if xrmin < lx[1] || xrmax > lx[2]
        error("Error in WKB: xrmin too small or xrmax too large!")
    end
    if yrmin < ly[1] || yrmax > ly[2]
        error("Error in WKB: yrmin too small or yrmax too large!")
    end
    if zrmin < lz[1] || zrmax > lz[2]
        error("Error in WKB: zrmin too small or zrmax too large!")
    end

    # Check if spectral-extent factors are set correctly.
    if sizex > 1 && fac_dk_init == 0.0
        error("Error in WKB: sizex > 1 && fac_dk_init == 0!")
    end
    if sizey > 1 && fac_dl_init == 0.0
        error("Error in WKB: sizey > 1 && fac_dl_init == 0!")
    end
    if sizez == 1 || fac_dm_init == 0.0
        error("Error in WKB: sizez == 1 || fac_dm_init == 0!")
    end

    # Set zonal index bounds and ray-volume count.
    if testcase == WKBMountainWave()
        ixmin = i0
        ixmax = i1
    else
        ixmin = max(i0, round(Int, (xr - lx[1] - dx / 2) / dx) + i0 - io)
        ixmax = min(i1, round(Int, (xr - lx[1] - dx / 2) / dx) + i0 - io)
    end
    if sizex == 1
        nxray = 1
    else
        nxray = nray_fac * nrxl * nrk_init
    end

    # Set meridional index bounds and ray-volume count.
    if testcase == WKBMountainWave()
        jymin = j0
        jymax = j1
    else
        jymin = max(j0, round(Int, (yr - ly[1] - dy / 2) / dy) + j0 - jo)
        jymax = min(j1, round(Int, (yr - ly[1] - dy / 2) / dy) + j0 - jo)
    end
    if sizey == 1
        nyray = 1
    else
        nyray = nray_fac * nryl * nrl_init
    end

    # Set vertical index bounds and ray-volume count.
    if testcase == WKBMountainWave()
        kzmin = k0 - 1
        kzmax = k0 - 1
        nzray = nray_fac * nrzl * nrm_init
    else
        kzmin = k0
        kzmax = k1
        nzray = nray_fac * nrzl * nrm_init
    end

    # Set maximum ray-volume count.
    nray_max = nxray * nyray * nzray * nwm

    # Set spectral dimension of ray-volume array.
    if nxray > 1
        nxray_wrk = 2 * nxray
    else
        nxray_wrk = 1
    end
    if nyray > 1
        nyray_wrk = 2 * nyray
    else
        nyray_wrk = 1
    end
    if nzray > 1
        nzray_wrk = 2 * nzray
    else
        nzray_wrk = 1
    end
    nray_wrk = nxray_wrk * nyray_wrk * nzray_wrk

    # Set number of surface ray volumes.
    n_sfc = nwm
    if nxray > 1
        n_sfc *= div(nxray, nray_fac)
    end
    if nyray > 1
        n_sfc *= div(nyray, nray_fac)
    end
    if nzray > 1
        n_sfc *= div(nzray, nray_fac)
    end

    # Initialize ray-volume arrays.
    nray = zeros(nxx, nyy, nzz)
    rays = Rays(nray_wrk, nxx, nyy, nzz)
    surface_indices = SurfaceIndices(n_sfc, nxx, nyy)
    increments = Increments(nray_wrk, nxx, nyy, nzz)
    integrals = Integrals(nxx, nyy, nzz)
    cgz_max = zeros(nxx, nyy, nzz)
    gwmomforce = Forces(nxx, nyy, nzz)

    # Initialize local arrays.
    omi_ini = zeros(nwm, nxx, nyy, nzz)
    wnk_ini = zeros(nwm, nxx, nyy, nzz)
    wnl_ini = zeros(nwm, nxx, nyy, nzz)
    wnm_ini = zeros(nwm, nxx, nyy, nzz)
    wad_ini = zeros(nwm, nxx, nyy, nzz)

    if wkb_mode == SteadyState() && testcase != WKBMountainWave()
        error(
            "Error in WKB: Steady state is implemented for WKBMountainWave only!",
        )
    end

    # TODO: correct size for zb?
    zb = zeros(nxx, nyy)

    if testcase == WKBMountainWave()
        # orographic_source!(
        #     namelists,
        #     constants,
        #     domain,
        #     grid,
        #     atmosphere,
        #     predictands,
        #     omi_ini,
        #     wnk_ini,
        #     wnl_ini,
        #     wnm_ini,
        #     wad_ini,
        #     zb,
        # )
    end

    # Initialize maximum horizontal group velocities.
    cgx_max = 0.0
    cgy_max = 0.0

    if ixmin <= ixmax && jymin <= jymax

        # Loop over all spatial cells with ray volumes.
        for kz in kzmin:kzmax, jy in jymin:jymax, ix in ixmin:ixmax
            iray = 0
            i_sfc = 0

            # Loop over all ray volumes within a spatial cell.
            for ix2 in 1:nrxl,
                ik in 1:nrk_init,
                jy2 in 1:nryl,
                jl in 1:nrl_init,
                kz2 in 1:nrzl,
                km in 1:nrm_init,
                iwm in 1:nwm

                # Set ray-volume indices.
                if testcase == WKBMountainWave()
                    i_sfc += 1

                    # Set surface indices.
                    surface_indices.ix2_sfc[i_sfc] = ix2
                    surface_indices.jy2_sfc[i_sfc] = jy2
                    surface_indices.kz2_sfc[i_sfc] = kz2
                    surface_indices.ik_sfc[i_sfc] = ik
                    surface_indices.jl_sfc[i_sfc] = jl
                    surface_indices.km_sfc[i_sfc] = km
                    surface_indices.iwm_sfc[i_sfc] = iwm

                    # Set surface ray-volume index.
                    if wad_ini[iwm, ix, jy, kz] == 0.0
                        surface_indices.ir_sfc[i_sfc, ix, jy] = -1
                        continue
                    else
                        iray += 1
                        surface_indices.ir_sfc[i_sfc, ix, jy] = iray
                    end
                else
                    iray += 1
                end

                # Set ray-volume positions.
                rays.x[iray, ix, jy, kz] =
                    (grid.x[io + ix] - 0.5 * dx + (ix2 - 0.5) * dx / nrxl)
                rays.y[iray, ix, jy, kz] =
                    (grid.y[jo + jy] - 0.5 * dy + (jy2 - 0.5) * dy / nryl)
                rays.z[iray, ix, jy, kz] = (
                    grid.ztfc[ix, jy, kz] - 0.5 * jac[ix, jy, kz] * dz +
                    (kz2 - 0.5) * jac[ix, jy, kz] * dz / nrzl
                )

                xr = rays.x[iray, ix, jy, kz]
                yr = rays.y[iray, ix, jy, kz]
                zr = rays.z[iray, ix, jy, kz]

                # Check if ray volume is too low.
                if zr < lz[1] - dz
                    error(
                        "Error in Rays: Ray volume",
                        iray,
                        "at",
                        ix,
                        jy,
                        kz,
                        "is too low!",
                    )
                end

                # Compute local stratification.
                n2r = stratification(zr, domain, grid, atmosphere, N2())

                # Set spatial extents.
                rays.dxray[iray, ix, jy, kz] = dx / nrxl
                rays.dyray[iray, ix, jy, kz] = dy / nryl
                rays.dzray[iray, ix, jy, kz] = jac[ix, jy, kz] * dz / nrzl

                wnk0 = wnk_ini[iwm, ix, jy, kz]
                wnl0 = wnl_ini[iwm, ix, jy, kz]
                wnm0 = wnm_ini[iwm, ix, jy, kz]

                # Ensure correct wavenumber extents.
                if testcase == WKBMountainWave() && sizex > 1
                    dk_ini_nd = fac_dk_init * sqrt(wnk0^2 + wnl0^2)
                end
                if testcase == WKBMountainWave() && sizey > 1
                    dl_ini_nd = fac_dl_init * sqrt(wnk0^2 + wnl0^2)
                end
                if wnm0 == 0.0
                    error("Error in WKB: wnm0 = 0!")
                else
                    dm_ini_nd = fac_dm_init * abs(wnm0)
                end

                # Set ray-volume wavenumbers.
                rays.k[iray, ix, jy, kz] =
                    (wnk0 - 0.5 * dk_ini_nd + (ik - 0.5) * dk_ini_nd / nrk_init)
                rays.l[iray, ix, jy, kz] =
                    (wnl0 - 0.5 * dl_ini_nd + (jl - 0.5) * dl_ini_nd / nrl_init)
                rays.m[iray, ix, jy, kz] =
                    (wnm0 - 0.5 * dm_ini_nd + (km - 0.5) * dm_ini_nd / nrm_init)

                # Set spectral extents.
                rays.dkray[iray, ix, jy, kz] = dk_ini_nd / nrk_init
                rays.dlray[iray, ix, jy, kz] = dl_ini_nd / nrl_init
                rays.dmray[iray, ix, jy, kz] = dm_ini_nd / nrm_init

                # Set phase-space volume.
                rays.area_xk[iray, ix, jy, kz] =
                    rays.dxray[iray, ix, jy, kz] * rays.dkray[iray, ix, jy, kz]
                rays.area_yl[iray, ix, jy, kz] =
                    rays.dyray[iray, ix, jy, kz] * rays.dlray[iray, ix, jy, kz]
                rays.area_zm[iray, ix, jy, kz] =
                    rays.dzray[iray, ix, jy, kz] * rays.dmray[iray, ix, jy, kz]
                pspvol = dm_ini_nd

                if sizex > 1
                    pspvol = pspvol * dk_ini_nd
                end
                if sizey > 1
                    pspvol = pspvol * dl_ini_nd
                end

                # Set phase-space wave-action density.
                if kz == sizez
                    rays.dens[iray, ix, jy, kz] = 0.0
                else
                    rays.dens[iray, ix, jy, kz] =
                        wad_ini[iwm, ix, jy, kz] / pspvol
                end

                # Set intrinsic frequency.
                rays.omega[iray, ix, jy, kz] = omi_ini[iwm, ix, jy, kz]

                # Interpolate winds to ray-volume position.
                uxr = meanflow(
                    xr,
                    yr,
                    zr,
                    namelists,
                    domain,
                    grid,
                    predictands,
                    U(),
                )
                vyr = meanflow(
                    xr,
                    yr,
                    zr,
                    namelists,
                    domain,
                    grid,
                    predictands,
                    V(),
                )
                wzr = meanflow(
                    xr,
                    yr,
                    zr,
                    namelists,
                    domain,
                    grid,
                    predictands,
                    W(),
                )

                wnrk = rays.k[iray, ix, jy, kz]
                wnrl = rays.l[iray, ix, jy, kz]
                wnrm = rays.m[iray, ix, jy, kz]
                wnrh = sqrt(wnrk^2 + wnrl^2)
                omir = rays.omega[iray, ix, jy, kz]

                # Compute maximum group velocities.
                cgirx = wnrk * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
                if abs(uxr + cgirx) > abs(cgx_max)
                    cgx_max = abs(uxr + cgirx)
                end
                cgiry = wnrl * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
                if abs(vyr + cgiry) > abs(cgy_max)
                    cgy_max = abs(vyr + cgiry)
                end
                cgirz =
                    -wnrm * (omir^2 - f_cor_nd^2) / (omir * (wnrh^2 + wnrm^2))
                if abs(wzr + cgirz) > abs(cgz_max)
                    cgz_max[ix, jy, kz] =
                        max(cgz_max[ix, jy, kz], abs[wzr + cgirz])
                end
            end

            # Set ray-volume count.
            nray[ix, jy, kz] = iray
            if iray > nray_wrk
                println(
                    "Error in Rays: nray[",
                    ix,
                    ", ",
                    jy,
                    ", ",
                    kz,
                    "] > nray_wrk =",
                    nray_wrk,
                )
                exit()
            end

            # Check if surface ray-volume count is correct.
            if testcase == WKBMountainWave()
                if i_sfc != n_sfc
                    println(
                        "Error in Rays: i_sfc =",
                        i_sfc,
                        "/= n_sfc =",
                        n_sfc,
                        "at [ix, jy, kz] = [",
                        ix,
                        ", ",
                        jy,
                        ", ",
                        kz,
                        "]",
                    )
                    exit()
                end
            end
        end
    end

    # Compute global ray-volume count.
    local_sum = sum(nray[i0:i1, j0:j1, (k0 - 1):k1])
    global_sum = MPI.Allreduce(local_sum, +, comm)

    # Print information.
    if master
        println("")
        println("Raytracer:")
        println("Global ray-volume count: ", global_sum)
        println("Maximum number of ray volumes per cell: ", nray_max)
    end

    return WKB(
        nxray,
        nyray,
        nzray,
        nxray_wrk,
        nyray_wrk,
        nzray_wrk,
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        increments,
        integrals,
        cgx_max,
        cgy_max,
        cgz_max,
        zb,
        gwmomforce,
    )
end
