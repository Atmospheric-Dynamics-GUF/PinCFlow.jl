struct IceNamelist{A <: AbstractIce}
    icesetup::A
end

function IceNamelist(; icesetup = NoIce())
    return IceNamelist(icesetup)
end
