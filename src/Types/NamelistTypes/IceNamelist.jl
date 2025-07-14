struct IceNamelist{A <: AbstractIce, B <: AbstractFloat}
     icesetup::A
     dt_ice::B
end

function IceNamelist(; icesetup = NoIce(), dt_ice = 1.0)
    return IceNamelist(icesetup, dt_ice)
end

