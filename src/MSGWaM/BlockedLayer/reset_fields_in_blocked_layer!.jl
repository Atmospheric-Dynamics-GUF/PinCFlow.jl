"""
```julia
reset_fields_in_blocked_layer!(state::State, choice::Symbol)::Nothing
```

Below the upper edge of the blocked layer, set all WKB fields of the kind specified by `choice` (either `:integrals` or `:tendencies`) to zero.

The fields in `state.wkb.integrals` are skipped, since they are taken care of in `include_blocked_flow_drag!`.

```julia
reset_fields_in_blocked_layer!(state::State, fields::AbstractWKBFields)::Nothing
```

Below the upper edge of the blocked layer, set all fields of `fields` to zero.

Each element `(i, j, k)` in the range `(i0:i1, j0:j1, k0:k1)` of each field is multiplied by `1 - blocked_fraction(state, i, j, k)`.

# Arguments

  - `state`: Model state.

  - `choice`: Chosen kind of WKB fields.

  - `fields`: Composite-type instance whose fields are to be modified.

# See also

  - [`PinCFlow.MSGWaM.BlockedLayer.blocked_fraction`](@ref)
"""
function reset_fields_in_blocked_layer! end

function reset_fields_in_blocked_layer!(state::State, choice::Symbol)::Nothing
    (; wkb, tracer, turbulence) = state
    (; tracer_setup) = state.namelists.tracer
    (; turbulence_scheme) = state.namelists.turbulence

    choice !== :tendencies && choice !== :integrals && error("Unknown choice!")

    fields = getfield(wkb, choice)
    tracer_fields = getfield(tracer, Symbol(:tracerwkb, choice))
    turbulence_fields = getfield(turbulence, Symbol(:turbulencewkb, choice))

    choice === :tendencies && reset_fields_in_blocked_layer!(state, fields)

    tracer_setup !== :NoTracer &&
        reset_fields_in_blocked_layer!(state, tracer_fields)

    turbulence_scheme !== :NoTurbulence &&
        reset_fields_in_blocked_layer!(state, turbulence_fields)

    nothing
end

function reset_fields_in_blocked_layer!(
    state::State,
    fields::AbstractWKBFields,
)::Nothing
    (; i0, i1, j0, j1, k0, k1) = state.domain

    for field_name in fieldnames(typeof(fields))
        field = getfield(fields, field_name)
        @fuse for k in k0:k1, j in j0:j1, i in i0:i1
            field[ii, jj, kk] *= (1 - blocked_fraction(state, i, j, k))
        end
    end

    nothing
end
