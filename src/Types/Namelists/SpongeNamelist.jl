struct SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractSponge,
    D <: Integer,
}
    spongelayer::A
    sponge_uv::A
    spongeheight::B
    spongealphaz_dim::B
    spongealphaz_fac::B
    unifiedsponge::A
    lateralsponge::A
    spongetype::C
    spongeorder::D
    cosmosteps::D
    relax_to_mean::A
    relaxation_period::B
    relaxation_amplitude::B
end

function SpongeNamelist(;
    spongelayer = false,
    sponge_uv = false,
    spongeheight = 5.0E-1,
    spongealphaz_dim = 1.0E-2,
    spongealphaz_fac = 1.0E+0,
    unifiedsponge = false,
    lateralsponge = false,
    spongetype = PolynomialSponge(),
    spongeorder = 1,
    cosmosteps = 1,
    relax_to_mean = true,
    relaxation_period = 0.0E+0,
    relaxation_amplitude = 0.0E+0,
)
    return SpongeNamelist(
        spongelayer,
        sponge_uv,
        spongeheight,
        spongealphaz_dim,
        spongealphaz_fac,
        unifiedsponge,
        lateralsponge,
        spongetype,
        spongeorder,
        cosmosteps,
        relax_to_mean,
        relaxation_period,
        relaxation_amplitude,
    )
end
