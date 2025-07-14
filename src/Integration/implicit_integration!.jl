function implicit_integration!(
    state::State,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    ntotalbicg::Integer,
    side::RHS,
    iout::Integer,
    machine_start_time::DateTime,
)
    (; master) = state.domain

    modify_compressible_wind!(state, *)

    set_boundaries!(state, BoundaryPredictands())

    save_backups!(state, :w)

    update!(state, dtstage, U(), side, Implicit(), 1.0)
    update!(state, dtstage, V(), side, Implicit(), 1.0)
    update!(state, dtstage, W(), side, Implicit(), 1.0)

    set_boundaries!(state, BoundaryPredictands())

    update!(state, dtstage, RhoP(), side, Implicit(), 1.0)

    set_boundaries!(state, BoundaryPredictands())

    (errflagbicg, niterbicg) = apply_corrector!(state, dtstage, 1.0, 1.0)

    if errflagbicg
        iout = write_output(state, time, iout, machine_start_time)
        if master
            println("Output last state into record", iout)
        end
        exit()
    end

    ntotalbicg += niterbicg

    modify_compressible_wind!(state, /)

    set_boundaries!(state, BoundaryPredictands())

    return ntotalbicg
end
