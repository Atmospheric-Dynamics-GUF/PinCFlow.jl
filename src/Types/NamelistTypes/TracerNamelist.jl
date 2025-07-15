struct TracerNamelist{A <: AbstractTracer, B <: Bool}
    tracersetup::A
    leading_order_impact::B
end

function TracerNamelist(;
    tracersetup = NoTracer(),
    leading_order_impact = true,
)
    return TracerNamelist(tracersetup, leading_order_impact)
end
