struct DiscretizationNamelist{
    A <: AbstractFloat,
    B <: Bool,
    C <: AbstractLimiter,
}
    cfl::A
    dtmin_dim::A
    dtmax_dim::A
    adaptive_time_step::B
    limitertype::C
end

function DiscretizationNamelist(;
    cfl = 5.0E-1,
    dtmin_dim = 1.0E-6,
    dtmax_dim = 1.0E+3,
    adaptive_time_step = true,
    limitertype = MCVariant(),
)
    return DiscretizationNamelist(
        cfl,
        dtmin_dim,
        dtmax_dim,
        adaptive_time_step,
        limitertype,
    )
end
