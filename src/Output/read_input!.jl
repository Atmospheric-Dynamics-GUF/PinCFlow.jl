"""
```julia
read_input!(state::State)
```

Read initial values for all prognostic variables from an HDF5 input file.

# Arguments

  - `state`: Model state.
"""
function read_input! end

function read_input!(state::State)
    (; ndx, ndy) = state.namelists.domain
    (; iin, input_file) = state.namelists.output
    (; model, test_case) = state.namelists.setting
    (; comm, nx, ny, nz, io, jo, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; lref, tref, rhoref, uref, thetaref) = state.constants
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands
    (; nray_max, nray, rays) = state.wkb
    (; rhostrattfc) = state.atmosphere

    # Determine dimensionality.
    dim = 1
    if ndx > 1
        dim += 1
    end
    if ndy > 1
        dim += 1
    end

    # Define slices.
    dk0 = ko == 0 ? 1 : 0
    (rr, ii, jj, kk, kkr) = (1:nray_max, i0:i1, j0:j1, k0:k1, (k0 - dk0):k1)
    (iid, jjd, kkd, kkrd) = (
        (io + 1):(io + nx),
        (jo + 1):(jo + ny),
        (ko + 1):(ko + nz),
        (ko + 2 - dk0):(ko + nz + 1),
    )

    # Open the file. Note: Fused in-place assignments cannot be used here!
    @ivy time = h5open(input_file, "r", comm) do file

        # Read the time.
        time = file["t"][iin] / tref

        # Read the density fluctuations.
        rhop[ii, jj, kk] = file["rhop"][iid, jjd, kkd, iin] ./ rhoref
        if model != Boussinesq()
            rho[ii, jj, kk] .= rhop[ii, jj, kk]
        end

        # Read the staggered zonal wind.
        u[ii, jj, kk] = file["us"][iid, jjd, kkd, iin] ./ uref

        # Read the staggered meridional wind.
        v[ii, jj, kk] = file["vs"][iid, jjd, kkd, iin] ./ uref

        # Read the staggered transformed vertical wind.
        w[ii, jj, kk] = file["wstfc"][iid, jjd, kkd, iin] ./ uref

        # Read the Exner-pressure fluctuations.
        pip[ii, jj, kk] = file["pip"][iid, jjd, kkd, iin]

        # Read the mass-weighted potential temperature.
        if model == Compressible()
            p[ii, jj, kk] = file["p"][iid, jjd, kkd, iin] ./ rhoref ./ thetaref
        end

        if !(typeof(state.namelists.tracer.tracer_setup) <: NoTracer)
            for field in fieldnames(TracerPredictands)
                getfield(state.tracer.tracerpredictands, field)[ii, jj, kk] =
                    file[string(field)][iid, jjd, kkd, iin] .*
                    (rhostrattfc[ii, jj, kk] .+ rho[ii, jj, kk])
            end
        end

        # Read ray-volume properties.
        if typeof(test_case) <: AbstractWKBTestCase
            for (output_name, field_name) in zip(
                ("xr", "yr", "zr", "dxr", "dyr", "dzr"),
                (:x, :y, :z, :dxray, :dyray, :dzray),
            )
                getfield(rays, field_name)[rr, ii, jj, kkr] =
                    file[output_name][rr, iid, jjd, kkrd, iin] ./ lref
            end

            for (output_name, field_name) in zip(
                ("kr", "lr", "mr", "dkr", "dlr", "dmr"),
                (:k, :l, :m, :dkray, :dlray, :dmray),
            )
                getfield(rays, field_name)[rr, ii, jj, kkr] =
                    file[output_name][rr, iid, jjd, kkrd, iin] .* lref
            end

            rays.dens[rr, ii, jj, kkr] =
                file["nr"][rr, iid, jjd, kkrd, iin] ./ rhoref ./ uref .^ 2 ./ tref ./ lref .^ dim

            # Determine nray.
            for k in kkr, j in jj, i in ii
                local_count = 0
                for r in rr
                    if rays.dens[r, i, j, k] == 0
                        continue
                    end
                    local_count += 1
                end
                nray[i, j, k] = local_count
            end
        end

        return time
    end

    return time
end
