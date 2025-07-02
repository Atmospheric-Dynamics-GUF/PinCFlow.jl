struct TurbulenceAuxiliaries{A <: AbstractFloat}
    tkebg::A
    ttebg::A
end

function TurbulenceAuxiliaries(constants::Constants)
    (; lref, tref) = constants

    # background tke and tte values
    tkebg = 0.1 * tref ^ 2.0 / lref ^ 2.0
    ttebg = 0.1 * tref ^ 2.0 / lref ^ 2.0

    return TurbulenceAuxiliaries(tkebg, ttebg)
end
