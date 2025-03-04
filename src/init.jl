function initialize_variables!(semi)
    (; equations, cache) = semi
    (; uRef) = equations
    (; var, rhoStrat, pStrat, thetaStrat, bvsStrat) = cache

    u = var.u
    v = var.v
    w = var.w
    rho = var.rho
    exner = var.exner
    rhop = var.rhop

    u0 = 10.0 # this should be a namelist parameter
    v0 = 0.0
    w0 = 0.0

    u .= u0 / uRef
    v .= v0 / uRef
    w .= w0 / uRef

    rho .= 0.0
    rhop .= 0.0
    exner .= 0.0
end

function initialize!(semi)
    initialize_atmosphere!(semi)

    initialize_variables!(semi)
end

function initialize_values(nx,
                           ny,
                           nz,
                           nbx,
                           nby,
                           nbz,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           zmin,
                           zmax;
                           model = "pseudo_incompressible",
                           maxIter = 5000,
                           tolcrit = "abs",
                           tolPoisson = 1e-8,
                           tolref = 1.0,)
    maxIter = 5000
    maxTime = 3.6e3
    tStepChoice = "cfl"
    dtMax_dim = 1.0
    cfl = 0.5e-1

    # Constants
    gamma = 1.4
    gamma_1 = gamma - 1.0
    kappa = (gamma - 1.0) / gamma
    kappaInv = 1.0 / kappa
    gammaInv = 1.0 / gamma

    # Reference quantities
    Rsp = 287.0
    g = 9.81
    rhoRef = 1.184
    pRef = 101325.0
    aRef = sqrt(pRef / rhoRef)
    uRef = aRef
    lRef = pRef / rhoRef / g
    tRef = lRef / aRef
    thetaRef = aRef^2.0 / Rsp
    Ma = uRef / aRef
    MaInv2 = 1.0 / Ma^2
    Fr = uRef / sqrt(g * lRef)
    kappa = (gamma - 1.0) / gamma
    sig = Ma^2 / Fr^2

    FRef = rhoRef * uRef^2 / lRef

    press0_dim = 1.0e5

    dt = dtMax_dim / tRef

    # isothermal atmosphere

    Temp0_dim = 300.0
    T0 = Temp0_dim / thetaRef
    N2 = Ma^2 / Fr^4 * kappa / T0
    NN = sqrt(N2)

    mu_viscous_dim = 0.0
    ReInv = mu_viscous_dim / (uRef * lRef)

    # Reynolds number
    ReInv = 0.0
    if false # TODO - Figure it out!
        ReInv = mu_viscous_dim / (uRef * lRef)
    end

    Re = if ReInv < 1.0e-20
        1.0e20
    else
        1.0 / ReInv
    end

    stretch_exponent = 1

    g_ndim = g / (uRef^2 / lRef)

    f_Coriolis_dim = 0.0

    grid = make_grid(nx = nx,
                     ny = ny,
                     nz = nz,
                     nbx = nbx,
                     nby = nby,
                     nbz = nbz,
                     xmin = xmin,
                     xmax = xmax,
                     ymin = ymin,
                     ymax = ymax,
                     zmin = zmin,
                     zmax = zmax,
                     lRef = lRef,
                     stretch_exponent = stretch_exponent)

    jac = Jacobian(grid)

    equations = (;
                 gamma,
                 gamma_1,
                 kappaInv,
                 gammaInv,
                 Rsp,
                 g,
                 rhoRef,
                 pRef,
                 aRef,
                 uRef,
                 lRef,
                 tRef,
                 thetaRef,
                 Ma,
                 Fr,
                 kappa,
                 sig,
                 press0_dim,
                 Temp0_dim,
                 T0,
                 N2,
                 NN,
                 mu_viscous_dim,
                 ReInv,
                 Re,
                 g_ndim,
                 MaInv2,
                 maxIter,
                 maxTime,
                 tStepChoice,
                 dtMax_dim,
                 cfl,
                 dt,
                 model,)

    # TODO: why not make Strat variables same size as var variables ??
    pStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                         (-nbx):(nx + nbx),
                         (-nby):(ny + nby),
                         -1:(nz + 2))

    thetaStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                             (-nbx):(nx + nbx),
                             (-nby):(ny + nby),
                             -1:(nz + 2))

    rhoStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                           (-nbx):(nx + nbx),
                           (-nby):(ny + nby),
                           -1:(nz + 2))

    bvsStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                           (-nbx):(nx + nbx),
                           (-nby):(ny + nby),
                           -1:(nz + 2))

    u = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz),
                    (-nbx):(nx + nbx),
                    (-nby):(ny + nby),
                    (-nbz):(nz + nbz))
    v, w, rho, exner, rhop = (copy(u), copy(u), copy(u), copy(u), copy(u))

    var = (; u, v, w, rho, exner, rhop)

    var0 = deepcopy(var)

    ndim = 3
    flux_u = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz, ndim),
                         (-nbx):(nx + nbx),
                         (-nby):(ny + nby),
                         (-nbz):(nz + nbz),
                         1:ndim)

    flux_u = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz, ndim),
                         (-nbx):(nx + nbx),
                         (-nby):(ny + nby),
                         (-nbz):(nz + nbz),
                         1:ndim)

    flux_v, flux_w, flux_rho, flux_exner, flux_rhop = (copy(flux_u), copy(flux_u),
                                                       copy(flux_u), copy(flux_u),
                                                       copy(flux_u))

    flux = (;
            u = flux_u,
            v = flux_v,
            w = flux_w,
            rho = flux_rho,
            exner = flux_exner,
            rhop = flux_rhop,)

    rhoTilde = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz, ndim,
                                 2),
                           (-nbx):(nx + nbx),
                           (-nby):(ny + nby),
                           (-nbz):(nz + nbz),
                           1:ndim,
                           0:1)

    uTilde, vTilde, wTilde, rhopTilde = (copy(rhoTilde), copy(rhoTilde), copy(rhoTilde),
                                         copy(rhoTilde))

    fluxDiffU = OffsetArray(zeros(2, 2), 0:1, 0:1)
    fluxDiffV = OffsetArray(zeros(2, 2), 0:1, 0:1)

    ## POISSON cache

    ac_b = Array{Float64, 3}(undef, nx, ny, nz)
    acv_b = Array{Float64, 3}(undef, nx, ny, nz)
    ach_b = Array{Float64, 3}(undef, nx, ny, nz)
    al_b = Array{Float64, 3}(undef, nx, ny, nz)
    ar_b = Array{Float64, 3}(undef, nx, ny, nz)
    ab_b = Array{Float64, 3}(undef, nx, ny, nz)
    af_b = Array{Float64, 3}(undef, nx, ny, nz)
    ad_b = Array{Float64, 3}(undef, nx, ny, nz)
    au_b = Array{Float64, 3}(undef, nx, ny, nz)
    aru_b = Array{Float64, 3}(undef, nx, ny, nz)
    ard_b = Array{Float64, 3}(undef, nx, ny, nz)
    alu_b = Array{Float64, 3}(undef, nx, ny, nz)
    ald_b = Array{Float64, 3}(undef, nx, ny, nz)
    afu_b = Array{Float64, 3}(undef, nx, ny, nz)
    afd_b = Array{Float64, 3}(undef, nx, ny, nz)
    abu_b = Array{Float64, 3}(undef, nx, ny, nz)
    abd_b = Array{Float64, 3}(undef, nx, ny, nz)
    auu_b = Array{Float64, 3}(undef, nx, ny, nz)
    add_b = Array{Float64, 3}(undef, nx, ny, nz)
    aruu_b = Array{Float64, 3}(undef, nx, ny, nz)
    ardd_b = Array{Float64, 3}(undef, nx, ny, nz)
    aluu_b = Array{Float64, 3}(undef, nx, ny, nz)
    aldd_b = Array{Float64, 3}(undef, nx, ny, nz)
    afuu_b = Array{Float64, 3}(undef, nx, ny, nz)
    afdd_b = Array{Float64, 3}(undef, nx, ny, nz)
    abuu_b = Array{Float64, 3}(undef, nx, ny, nz)
    abdd_b = Array{Float64, 3}(undef, nx, ny, nz)

    matVec = Array{Float64, 3}(undef, nx, ny, nz)

    v_pc = Array{Float64, 3}(undef, nx, ny, nz)

    dp_ = zeros(nx + 2, ny + 2, nz + 2)
    dp = OffsetArray(dp_, 0:(nx + 1), 0:(ny + 1), 0:(nz + 1))

    s_pc = zeros(nx, ny, nz)
    q_pc = zeros(nx, ny, nz)
    p_pc = zeros(nx, ny)

    kr_sp_tfc_ = zeros(nx + 2, ny + 2, nz + 2)
    kr_sp_tfc = OffsetArray(kr_sp_tfc_, 0:(nx + 1), 0:(ny + 1), 0:(nz + 1))

    kr_sp_w_tfc = copy(kr_sp_tfc)

    rhs_bicg = zeros(nx, ny, nz)

    s_aux_field_lin_opr_ = zeros(nx + 2, ny + 2, nz + 2)
    s_aux_field_lin_opr = OffsetArray(s_aux_field_lin_opr_, 0:(nx + 1), 0:(ny + 1),
                                      0:(nz + 1))

    sum_local_bicg = zeros(nz)
    sum_global_bicg = zeros(nz)

    sol_bicg = zeros(nx, ny, nz)

    r_vm = zeros(nx, ny)

    p_bicg = zeros(nx, ny, nz)
    r0_bicg = zeros(nx, ny, nz)
    rOld_bicg = zeros(nx, ny, nz)
    r_bicg = zeros(nx, ny, nz)
    s_bicg = zeros(nx, ny, nz)
    b_bicg = zeros(nx, ny, nz)
    t_bicg = zeros(nx, ny, nz)
    v_bicg = zeros(nx, ny, nz)
    matVec_bicg = zeros(nx, ny, nz)
    v_pc_bicg = zeros(nx, ny, nz)

    corX = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                       (-nbx):(nx + nbx),
                       (-nby):(ny + nby),
                       -1:(nz + 2))
    corY = copy(corX)

    cache = (;
             var,
             var0,
             flux,
             uTilde,
             vTilde,
             wTilde,
             rhoTilde,
             rhopTilde,
             pStrat,
             thetaStrat,
             rhoStrat,
             bvsStrat,
             jac,
             fluxDiffU,
             fluxDiffV,
             ac_b,
             acv_b,
             ach_b,
             al_b,
             ar_b,
             ab_b,
             af_b,
             ad_b,
             au_b,
             aru_b,
             ard_b,
             alu_b,
             ald_b,
             afu_b,
             afd_b,
             abu_b,
             abd_b,
             auu_b,
             add_b,
             aruu_b,
             ardd_b,
             aluu_b,
             aldd_b,
             afuu_b,
             afdd_b,
             abuu_b,
             abdd_b,
             matVec,
             dp,
             v_pc,
             s_pc,
             q_pc,
             p_pc,
             kr_sp_tfc,
             kr_sp_w_tfc,
             rhs_bicg,
             s_aux_field_lin_opr,
             sum_local_bicg,
             sum_global_bicg,
             sol_bicg,
             r_vm,
             p_bicg,
             r0_bicg,
             rOld_bicg,
             r_bicg,
             s_bicg,
             b_bicg,
             t_bicg,
             v_bicg,
             matVec_bicg,
             v_pc_bicg,
             corX,
             corY,)
    #  f_cor_nd =  OffsetArray(zeros(Float64, 1:ny), -1:ny-1)
    f_cor_nd = 0.0
    boundary_x = boundary_y = PeriodicBC()
    boundary_z = SolidWallBC()
    boundary_conditions = (; boundary_x, boundary_y, boundary_z)

    parameters = (; maxIter, tolcrit, tolPoisson, tolref)

    met = MetricTensor(cache, grid, equations)

    sponge = SinusodialSponge(grid)

    sponge_layer = true

    sponge_uv = false

    return (;
            grid,
            equations,
            cache,
            boundary_conditions,
            parameters,
            met,
            sponge,
            sponge_layer,
            sponge_uv,)
end

function time_discretization(semi, dt)
    (; cache) = semi
    (; var, flux) = cache

    betaRK = (1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0)

    alphaRK = (0.0, -5.0 / 9.0, -153.0 / 128.0)

    stepFrac = (1 / 3, 5 / 12, 1 / 4)

    dRho = copy(var.u)
    dRhop = copy(var.u)
    dMom = copy(flux.u)
    usave = copy(var.u)
    vsave = copy(var.v)
    wsave = copy(var.w)
    uOld = copy(var.u)
    vOld = copy(var.v)
    wOld = copy(var.w)
    rhoOld = copy(var.u)
    rhopOld = copy(var.u)
    flux0 = deepcopy(flux)
    var0 = deepcopy(var)
    var1 = deepcopy(var)
    return (;
            alphaRK,
            betaRK,
            dRho,
            dRhop,
            dMom,
            usave,
            vsave,
            wsave,
            uOld,
            vOld,
            wOld,
            rhoOld,
            rhopOld,
            dt,
            flux0,
            var0,
            var1,
            stepFrac,)
end
