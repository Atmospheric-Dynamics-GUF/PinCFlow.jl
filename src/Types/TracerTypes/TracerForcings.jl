struct TracerForcings{A <: TracerGWImpact}
    chiq0::A
end

function TracerForcings(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer

    return TracerForcings(namelists, domain, tracersetup)
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracersetup::NoTracer,
)
    return TracerForcings(TracerGWImpact(0, 0, 0))
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracersetup::AbstractTracer,
)
    (; testcase) = namelists.setting

    return TracerForcings(namelists, domain, testcase)
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    testcase::AbstractTestCase,
)
    return TracerForcings(TracerGWImpact(0, 0, 0))
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    testcase::AbstractWKBTestCase,
)
    (; nxx, nyy, nzz) = domain

    (; model) = namelists.setting 
    (; leading_order_impact) = namelists.tracer

    if leading_order_impact && model == PseudoIncompressible()
        error("Cannot calculate leading_order_impact in TFC with PI model")
    end

    return TracerForcings(TracerGWImpact(nxx, nyy, nzz))
end