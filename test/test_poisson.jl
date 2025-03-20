using PinCFlow

include("util.jl")
function tensor_norms(tensor)
  return get_norms.((
    tensor.ac_b,
    tensor.acv_b,
    tensor.ach_b,
    tensor.al_b,
    tensor.ar_b,
    tensor.ab_b,
    tensor.af_b,
    tensor.ad_b,
    tensor.au_b,
    tensor.aru_b,
    tensor.ard_b,
    tensor.alu_b,
    tensor.ald_b,
    tensor.afu_b,
    tensor.afd_b,
    tensor.abu_b,
    tensor.abd_b,
    tensor.auu_b,
    tensor.add_b,
    tensor.aruu_b,
    tensor.ardd_b,
    tensor.aluu_b,
    tensor.aldd_b,
    tensor.afuu_b,
    tensor.afdd_b,
    tensor.abuu_b,
  ))
end

domain = PinCFlow.DomainNamelist(; sizex = 300, sizey = 1, sizez = 100)
namelists = PinCFlow.Namelists(; domain = domain)
state = PinCFlow.State(namelists)

normalizer = 0.1
function set_lin_array_normalizer!(arr)
  return set_lin_array!(arr; normalizer = normalizer)
end

set_lin_array_normalizer!.((
  state.atmosphere.pstrattfc,
  state.atmosphere.rhostrattfc,
  state.atmosphere.bvsstrattfc,
  state.sponge.kr_sp_tfc,
  state.sponge.kr_sp_w_tfc,
  state.variables.predictands.rho,
))

set_lin_array_normalizer!.((
  state.grid.topography_surface,
  state.grid.ztildes,
  state.grid.zs,
))

PinCFlow.compute_operator!(state, 1.0, PinCFlow.EXPL(), 1.0)

norms = tensor_norms(state.poisson.tensor)

@test norms[1] ==
      (2.3241474438706848e11, 1.3419005488701441e9, 7.788720155394001e6)
@test norms[2] == (2.6939016189804436e10, 1.559938266764099e8, 939522.101139049)
@test norms[3] ==
      (2.0547736256932065e11, 1.1863241176531262e9, 6.960914045460912e6)
@test norms[4] ==
      (1.0273961573006381e11, 5.931685732645231e8, 3.517691919706872e6)
@test norms[5] ==
      (1.0273546378307071e11, 5.931446233714064e8, 3.519857949212931e6)
@test norms[12] == (9.42463084628637e8, 7.26926739962243e6, 234554.58480521778)

# @testset "preCond expl" begin
#   semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
#   (; cache, grid) = semi
#   (;
#     p_bicg,
#     p_pc,
#     v_pc,
#     au_b,
#     ad_b,
#     ac_b,
#     al_b,
#     ar_b,
#     af_b,
#     ab_b,
#     ach_b,
#     acv_b,
#     aru_b,
#     ard_b,
#     alu_b,
#     ald_b,
#     afu_b,
#     afd_b,
#     abu_b,
#     abd_b,
#     auu_b,
#     add_b,
#     aruu_b,
#     ardd_b,
#     aluu_b,
#     aldd_b,
#     afuu_b,
#     afdd_b,
#     abuu_b,
#     abdd_b,
#     q_pc,
#     p_pc,
#   ) = cache

#   normalizer = 0.01 # To get smaller values where arithmetic is more accurate
#   set_lin_array_normalizer!(arr) = set_lin_array!(arr; normalizer = normalizer)
#   set_lin_array_normalizer!.((
#     p_bicg,
#     p_pc,
#     v_pc,
#     au_b,
#     ad_b,
#     ac_b,
#     al_b,
#     ar_b,
#     af_b,
#     ab_b,
#     ach_b,
#     acv_b,
#     aru_b,
#     ard_b,
#     alu_b,
#     ald_b,
#     afu_b,
#     afd_b,
#     abu_b,
#     abd_b,
#     auu_b,
#     add_b,
#     aruu_b,
#     ardd_b,
#     aluu_b,
#     aldd_b,
#     afuu_b,
#     afdd_b,
#     abuu_b,
#     abdd_b,
#     q_pc,
#     p_pc,
#   ))

#   PinCFlow.preCond(p_bicg, v_pc, "expl", semi)
#   ref_filename = "$(pincflow_test_dir())/poisson_fortran_data/preCond/sOut.txt.gz"
#   test_file(ref_filename, v_pc; tol_l1 = 1e-16, tol_linf = 1e-16)
# end

# @testset "calc_RHS and poissonSolver" begin
#   semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
#   initialize_atmosphere!(semi)
#   initialize_variables!(semi)
#   (; cache, grid) = semi
#   rhs = cache.rhs_bicg
#   (; pStrat, rhoStrat, var) = cache
#   (; topography_surface) = grid

#   normalizer = 0.01 # To get smaller values where arithmetic is more accurate
#   set_lin_array_normalizer!(arr) = set_lin_array!(arr; normalizer = normalizer)
#   set_lin_array_normalizer!.((pStrat, rhoStrat, topography_surface))

#   set_lin_array_normalizer!(var.u)
#   set_lin_array_normalizer!(var.v)
#   set_lin_array_normalizer!(var.w)

#   PinCFlow.calc_RHS(rhs, semi, 1.0)
#   ref_filename = "$(pincflow_test_dir())/poisson_fortran_data/calc_RHS/rhs.txt.gz"
#   test_file(ref_filename, rhs; tol_l1 = 2e-13, tol_linf = 2e-13)
# end

# @testset "Poisson solver" begin
#   semi = initialize_values(30, 1, 10, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
#   initialize_atmosphere!(semi)
#   initialize_variables!(semi)
#   (; cache, grid) = semi
#   rhs = cache.rhs_bicg
#   (; pStrat, rhoStrat, var, kr_sp_tfc, kr_sp_w_tfc) = cache
#   (; topography_surface) = grid

#   normalizer = 0.01 # To get smaller values where arithmetic is more accurate
#   set_lin_array_normalizer!(arr) = set_lin_array!(arr; normalizer = normalizer)
#   set_lin_array_normalizer!.((pStrat, rhoStrat, topography_surface))

#   set_lin_array_normalizer!(var.u)
#   set_lin_array_normalizer!(var.v)
#   set_lin_array_normalizer!(var.w)

#   PinCFlow.calc_RHS(rhs, semi, 1.0)
#   dt = 0.3
#   errFlagBicg = false
#   nIter = 0
#   facprs = 1.0
#   facray = 1.0
#   PinCFlow.poissonSolver(
#     rhs,
#     semi,
#     dt,
#     errFlagBicg,
#     nIter,
#     "expl",
#     facray,
#     facprs,
#   )

#   test_arr(
#     cache.dp,
#     6954.114011644017,
#     406.0728580408809,
#     32.1315910002854;
#     tol = 1e-9,
#   )

#   PinCFlow.poissonSolver(
#     rhs,
#     semi,
#     dt,
#     errFlagBicg,
#     nIter,
#     "impl",
#     facray,
#     facprs,
#   )

#   test_arr(
#     cache.dp,
#     6949.851010678487,
#     405.8250898219831,
#     32.1214829208721;
#     tol = 1e-9,
#   )

#   PinCFlow.pressureBoundaryCondition(semi)

#   test_arr(
#     cache.dp,
#     26794.367327608583,
#     799.5740338430224,
#     32.1214829208721;
#     tol = 1e-9,
#   )

#   PinCFlow.correctorStep(semi, dt, "expl", facray, facprs)

#   test_arr(
#     var.exner,
#     26794.367327608583,
#     799.5740338430224,
#     32.1214829208721;
#     tol = 1e-9,
#   )
#   test_arr(
#     var.u,
#     2996.1961128135836,
#     250.73352464698718,
#     64.47576033368075;
#     tol = 1e-9,
#   )

#   PinCFlow.correctorStep(semi, dt, "impl", facray, facprs)

#   set_lin_array_normalizer!.((kr_sp_tfc, kr_sp_w_tfc))
#   test_arr(
#     var.exner,
#     53588.734655217166,
#     1599.148067686045,
#     64.2429658417442;
#     tol = 1e-9,
#   )
#   test_arr(kr_sp_tfc, 368.64, 11.368166079012056, 0.54; tol = 1e-9)
#   test_arr(kr_sp_w_tfc, 368.64, 11.368166079012056, 0.54; tol = 1e-9)
#   test_arr(var.u, 1534.76, 23.24188460516906, 0.6; tol = 1e-9)
#   test_arr(
#     cache.corX,
#     1643.1243369208187,
#     248.88860795404267,
#     64.05576033368075;
#     tol = 1e-9,
#   )
#   test_arr(
#     cache.corY,
#     0.15327530873362505,
#     0.013421375877617665,
#     0.0033422640397776215;
#     tol = 1e-9,
#   )
#   test_arr(
#     var.exner,
#     53588.734655217166,
#     1599.148067686045,
#     64.2429658417442;
#     tol = 1e-9,
#   )

#   PinCFlow.bicgstab(rhs, dt, semi, cache.sol_bicg, nIter, errFlagBicg, "expl")

#   test_arr(
#     cache.sol_bicg,
#     1142.0651583477945,
#     66.0122418332421,
#     4.178712816389702;
#     tol = 1e-9,
#   )
#   test_arr(
#     cache.p_bicg,
#     5.676571678854043e-6,
#     3.743420434190437e-7,
#     3.6087341096433547e-8;
#     tol = 1e-9,
#   )

#   PinCFlow.bicgstab(rhs, dt, semi, cache.sol_bicg, nIter, errFlagBicg, "impl")

#   test_arr(
#     cache.sol_bicg,
#     1142.0651583477945,
#     66.0122418332421,
#     4.178712816389702;
#     tol = 1e-9,
#   )
#   test_arr(
#     cache.p_bicg,
#     5.676571678854043e-6,
#     3.743420434190437e-7,
#     3.6087341096433547e-8;
#     tol = 1e-9,
#   )
# end
