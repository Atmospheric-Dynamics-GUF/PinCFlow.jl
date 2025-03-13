abstract type AbstractSponge end
struct ExponentialSponge <: AbstractSponge end
struct COSMOSponge <: AbstractSponge end
struct PolynomialSponge <: AbstractSponge end
struct SinusoidalSponge <: AbstractSponge end

struct SpongeNamelist{
    A<:Bool,
    B<:AbstractFloat,
    C<:AbstractSponge,
    D<:Integer,
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
  spongeheight = 0.33,
  spongealphaz_dim = 0.01,
  spongealphaz_fac = 0.01,
  unifiedsponge = false,
  lateralsponge = false,
  spongetype = PolynomialSponge(),
  spongeorder = 1,
  cosmosteps = 1,
  relax_to_mean = true,
  relaxation_period = 0.0,
  relaxation_amplitude = 0.0,
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
