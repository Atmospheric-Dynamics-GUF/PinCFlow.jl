struct TracerNamelist{A <: AbstractTracer}
    tracersetup::A
end

function TracerNamelist(; tracersetup = NoTracer())
    return TracerNamelist(tracersetup)
end
