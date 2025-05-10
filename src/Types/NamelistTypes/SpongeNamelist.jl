struct SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractSponge,
    D <: Integer,
    E <: NTuple{3, <:AbstractFloat},
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
    perturbation_period::B
    perturbation_amplitude::B
    relaxation_wind::E
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
    perturbation_period = 0.0E+0,
    perturbation_amplitude = 0.0E+0,
    relaxation_wind = (0.0E+0, 0.0E+0, 0.0E+0),
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
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    )
end
