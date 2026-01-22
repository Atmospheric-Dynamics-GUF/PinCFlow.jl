"""
```julia
IceNamelist{A <: AbstractIce}
```

Namelist for the inclusion of ice physics.

```julia
IceNamelist(; ice_setup::AbstractIce = NoIce())::IceNamelist
```

Construct an `IceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `ice_setup::A`: General ice-physics configuration.
"""
struct IceNamelist{A <: AbstractIce, B <: AbstractFloat, C <: Integer, D <: AbstractCloudCover, E <: Bool, F <: AbstractIceTestCase}
<<<<<<< HEAD
     tau_q_sink::B # time scale for q sink (s)
     icesetup::A
=======
     ice_setup::A
>>>>>>> origin/ice
     dt_ice::B
     nscx :: C
     nscy :: C
     nscz :: C
     cloudcover :: D
     parameterized_nucleation :: E
     parameterized_sgs_q :: E
     constant_advection :: E
     hor_adv_vel :: NTuple{2, <:AbstractFloat}
     ice_test_case :: F
end

function IceNamelist(; tau_q_sink = 3.0e-11, ice_setup = NoIce(), dt_ice = 1.0, nscx = 1, nscy = 1, nscz = 1, cloudcover= CloudCoverOff(), parameterized_nucleation = false, parameterized_sgs_q = false, constant_advection = false, hor_adv_vel = (0.0, 0.0), ice_test_case = NoIceTestCase())
    return IceNamelist(tau_q_sink, ice_setup, dt_ice, nscx, nscy, nscz, cloudcover, parameterized_nucleation, parameterized_sgs_q, constant_advection, hor_adv_vel, ice_test_case)
end

