
struct DomainParameters{I<:Integer,VOF<:Vector{<:AbstractFloat}}
    sizex::I
    sizey::I
    sizez::I
    nbx::I
    nby::I
    nbz::I
    lx_dim::VOF
    ly_dim::VOF
    lz_dim::VOF
    nprocx::I
    nprocy::I
end

function DomainParameters(;
    sizex=150,
    sizey=1,
    sizez=50,
    nbx=3,
    nby=3,
    nbz=3,
    lx_dim=[0.0E0, 6.0E+4],
    ly_dim=[0.0E0, 4.0E+4],
    lz_dim=[0.0E0, 2.0E+4],
    nprocx=1,
    nprocy=1)
    DomainParameters(sizex, sizey, sizez, nbx, nby, nbz,
        lx_dim, ly_dim, lz_dim, nprocx,
        nprocy)
end

struct OutputParameters{S<:AbstractString,I<:Integer,F<:AbstractFloat,
    VOS<:Vector{S}}
    atmvarout::VOS
    prepare_restart::Bool
    restart::Bool
    iin::I
    runname::S
    outputtype::S
    noutput::I
    maxiter::I
    outputtimediff::F
    maxtime::F
    fancy_namelists::Bool
end

function OutputParameters(;
    atmvarout=["w"],
    prepare_restart=false,
    restart=false,
    iin=-1,
    runname="mountainwave",
    outputtype="time",
    noutput=1,
    maxiter=1,
    outputtimediff=3.6E3,
    maxtime=3.6E3,
    fancy_namelist=true)
    OutputParameters(atmvarout, prepare_restart, restart, iin, runname, outputtype, noutput,
        maxiter, outputtimediff, maxtime, fancy_namelist)
end

struct DebugParameters{F<:AbstractFloat}
    dtmin_dim::F
end

DebugParameters(; dtmin_dim=1.0E-5) = DebugParameters(dtmin_dim)

struct TestCaseParameters
    testcase::AbstractString
end
TestCaseParameters(; testcase="mountainwave") = TestCaseParameters(testcase)

struct ModelParameters
    model::AbstractString
end
ModelParameters(; model="pseudo_incompressible") = ModelParameters(model)

struct DiscretizationParameters{F<:AbstractFloat,S<:AbstractString}
    cfl::F
    dtmax_dim::F
    tstepchoice::S
    limitertype1::S
end

function DiscretizationParameters(;
    cfl=5.0E-1,
    dtmax_dim=6.0E+1,
    tstepchoice="cfl",
    limitertype1="MCVariant")
    DiscretizationParameters(cfl, dtmax_dim, tstepchoice, limitertype1)
end

struct PoissonSolverParameters{F<:AbstractFloat,I<:Integer,S<:AbstractString}
    tolpoisson::F
    maxiter::I
    preconditioner::S
    dtau::F
    maxiteradi::I
    initialcleaning::Bool
    correctmomentum::Bool
    tolcrit::S
    tolref::F
end

function PoissonSolverParameters(;
    tolpoisson=1.0E-8, maxiterpoisson=5000,
    preconditioner="yes", dtau=4.0E0, maxiteradi=2,
    initalcleaning=true, correctmomentum=true,
    tolcrit="abs", tolref=1.0)

    PoissonSolverParameters(tolpoisson, maxiterpoisson, preconditioner, dtau, maxiteradi,
        initalcleaning, correctmomentum, tolcrit, tolref)
end

struct AtmosphereParameters{F<:AbstractFloat,S<:AbstractString,VOF<:Vector{F}}
    specifyreynolds::Bool
    reinv::F
    mu_viscous_dim::F
    background::S
    temp0_dim::F
    press0_dim::F
    backgroundflow_dim::VOF
    f_coriolis_dim::F
    corset::S
end

function AtmosphereParameters(;
    specifyreynolds=false,
    reinv=0.0,
    mu_viscous_dim=0.0,
    background="isothermal",
    temp0_dim=3.0E2,
    press0_dim=1.0E5,
    backgroundflow_dim=[1.0E1, 0.0, 0.0],
    f_coriolis_dim=0.0,
    corset="constant")
    AtmosphereParameters(specifyreynolds, reinv, mu_viscous_dim, background, temp0_dim,
        press0_dim, backgroundflow_dim, f_coriolis_dim, corset)
end

struct TopographyParameters{F<:AbstractFloat,I<:Integer}
    mountainheight_dim::F
    mountainwidth_dim::F
    mountain_case::I
    range_factor::F
    spectral_modes::I
    envelope_reduction::F
    stretch_exponent::F
end

function TopographyParameters(;
    mountainheight_dim=4.0E2,
    mountainwidth_dim=1.0E3,
    mountain_case=3,
    range_factor=10.0,
    spectral_modes=1,
    envelope_reduction=0.0,
    stretch_exponent=1.0,)
    TopographyParameters(mountainheight_dim, mountainwidth_dim, mountain_case, range_factor,
        spectral_modes, envelope_reduction, stretch_exponent)
end

struct BoundaryParameters{F<:AbstractFloat,S<:AbstractString,I<:Integer}
    spongelayer::Bool
    sponge_uv::Bool
    spongeheight::F
    spongealphaz_dim::F
    spongealphaz_fac::F
    unifiedsponge::Bool
    lateralsponge::Bool
    spongetype::S
    spongeorder::I
    cosmosteps::I
    relax_to_mean::Bool
    relaxation_period::F
    relaxation_amplitude::F
    xboundary::S
    yboundary::S
    zboundary::S
end

function BoundaryParameters(;
    spongelayer=true,
    sponge_uv=false,
    spongeheight=5.0E-1,
    spongealphaz_dim=1.0E-2,
    spongealphaz_fac=1.0,
    unifiedsponge=true,
    lateralsponge=true,
    spongetype="sinusodial",
    spongeorder=1,
    cosmosteps=1,
    relax_to_mean=true,
    relaxation_period=0.0,
    relaxation_amplitude=0.0,
    xboundary="periodic",
    yboundary="periodic",
    zboundary="solid_wall",)
    BoundaryParameters(spongelayer, sponge_uv, spongeheight, spongealphaz_dim,
        spongealphaz_fac, unifiedsponge, lateralsponge, spongetype,
        spongeorder, cosmosteps, relax_to_mean, relaxation_period,
        relaxation_amplitude, xboundary, yboundary, zboundary)
end

struct Parameters
    domain::DomainParameters
    output::OutputParameters
    debug::DebugParameters
    testcase::TestCaseParameters
    model::ModelParameters
    discretization::DiscretizationParameters
    poisson::PoissonSolverParameters
    atmosphere::AtmosphereParameters
    topography::TopographyParameters
    boundaries::BoundaryParameters
end

function default_parameters()
    return Parameters(
        DomainParameters(),
        OutputParameters(),
        DebugParameters(),
        TestCaseParameters(),
        ModelParameters(),
        DiscretizationParameters(),
        PoissonSolverParameters(),
        AtmosphereParameters(),
        TopographyParameters(),
        BoundaryParameters()
    )
end
