using PinCFlow_dev
using GZip

function set_lin_array!(arr::AbstractArray; normalizer::Real = 1.0)
    for I in CartesianIndices(arr)
        I_tuple = get_ijk(I) # I_tuple = (i, j, k)
        arr[I_tuple...] = (sum(I_tuple) + 10) * normalizer # arr[i, j, k] = (i + j + k)^2 + 10
    end
end

get_ijk(I::CartesianIndex{3}) = I[1], I[2], I[3]
get_ijk(I::CartesianIndex{2}) = I[1], I[2]
get_ijk(I::CartesianIndex{1}) = I[1]

macro get_var_name(var)
    return QuoteNode(var)  # Returns the variable name as a Symbol
end

function test_file(file, data; tol_l1 = 1e-12, tol_linf = 1e-10)
    error_l1 = 0.0
    error_linf = 0.0
    GZip.open(file) do io
        for k in axes(data, 3)
            for j in axes(data, 2)
                for i in axes(data, 1)
                    s = readline(io)
                    ref_val = parse(Float64, s)
                    err = abs(data[i, j, k] - ref_val)
                    error_linf = max(error_linf, err)
                    error_l1 += err
                end
            end
        end
    end
    error_l1 /= length(data)
    @test error_l1 < tol_l1
    @test error_linf < tol_linf
end

## Tests begin

@testset "val_PsIn expl" begin
    semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)

    (; equations, cache, grid, met) = semi
    (; nx, ny, nz, dx, dy, dz) = grid
    (; pStrat, rhoStrat, jac, var, kr_sp_tfc, kr_sp_w_tfc, bvsStrat) = cache

    # Set "random values"

    set_lin_array!.((pStrat, rhoStrat, kr_sp_tfc, kr_sp_w_tfc, bvsStrat, cache.var.rho))

    set_lin_array!.((grid.topography_surface, grid.zTildeS, grid.zS))

    # Call function to be tested
    PinCFlow_dev.val_PsIn(semi, 1.0, "expl", 1.0)

    # Arrange values to be tested
    (; ac_b, acv_b, ach_b, al_b, ar_b, ab_b, af_b, ad_b, au_b, aru_b, ard_b,
    alu_b, ald_b, afu_b, afd_b, abu_b, abd_b, auu_b, add_b, aruu_b, ardd_b,
    aluu_b, aldd_b, afuu_b, afdd_b, abuu_b, abdd_b,) = cache

    test_cache = (;
                  ac_b,
                  acv_b, ach_b, al_b, ar_b, ab_b, af_b, ad_b, au_b, aru_b,
                  ard_b, alu_b, ald_b, afu_b, afd_b, abu_b, abd_b, auu_b,
                  add_b, aruu_b, ardd_b, aluu_b, aldd_b, afuu_b, afdd_b,
                  abuu_b,)

    for field in fieldnames(typeof(test_cache))
        val = getfield(test_cache, field)
        true_arr = getfield(cache, field)
        ref_filename = "$(pincflow_test_dir())/poisson_fortran_data/val_PsIn/$(field).txt.gz"
        test_file(ref_filename, val)
    end
end

@testset "preCond expl" begin
    semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    (; cache, grid) = semi
    (; p_bicg, p_pc, v_pc, au_b, ad_b, ac_b, al_b, ar_b, af_b, ab_b, ach_b,
    acv_b, aru_b, ard_b, alu_b, ald_b, afu_b, afd_b, abu_b, abd_b, auu_b,
    add_b, aruu_b, ardd_b, aluu_b, aldd_b, afuu_b, afdd_b, abuu_b, abdd_b,
    q_pc, p_pc,) = cache

    normalizer = 0.01 # To get smaller values where arithmetic is more accurate
    set_lin_array_normalizer!(arr) = set_lin_array!(arr, normalizer = normalizer)
    set_lin_array_normalizer!.((p_bicg, p_pc, v_pc, au_b, ad_b, ac_b, al_b, ar_b, af_b,
                                ab_b, ach_b,
                                acv_b, aru_b, ard_b, alu_b, ald_b, afu_b, afd_b, abu_b,
                                abd_b, auu_b,
                                add_b, aruu_b, ardd_b, aluu_b, aldd_b, afuu_b, afdd_b,
                                abuu_b, abdd_b, q_pc,
                                p_pc))

    PinCFlow_dev.preCond(p_bicg, v_pc, "expl", semi)
    ref_filename = "$(pincflow_test_dir())/poisson_fortran_data/preCond/sOut.txt.gz"
    test_file(ref_filename, v_pc, tol_l1 = 1e-16, tol_linf = 1e-16)
end

@testset "calc_RHS and poissonSolver" begin
    semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    initialize_atmosphere!(semi)
    initialize_variables!(semi)
    (; cache, grid) = semi
    rhs = cache.rhs_bicg
    (; pStrat, rhoStrat, var) = cache
    (; topography_surface) = grid

    normalizer = 0.01 # To get smaller values where arithmetic is more accurate
    set_lin_array_normalizer!(arr) = set_lin_array!(arr, normalizer = normalizer)
    set_lin_array_normalizer!.((pStrat, rhoStrat, topography_surface))

    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)

    PinCFlow_dev.calc_RHS(rhs, semi, 1.0)
    ref_filename = "$(pincflow_test_dir())/poisson_fortran_data/calc_RHS/rhs.txt.gz"
    test_file(ref_filename, rhs, tol_l1 = 2e-13, tol_linf = 2e-13)
end

@testset "Poisson solver" begin
    semi = initialize_values(30, 1, 10, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    initialize_atmosphere!(semi)
    initialize_variables!(semi)
    (; cache, grid) = semi
    rhs = cache.rhs_bicg
    (; pStrat, rhoStrat, var, kr_sp_tfc, kr_sp_w_tfc) = cache
    (; topography_surface) = grid

    normalizer = 0.01 # To get smaller values where arithmetic is more accurate
    set_lin_array_normalizer!(arr) = set_lin_array!(arr, normalizer = normalizer)
    set_lin_array_normalizer!.((pStrat, rhoStrat, topography_surface))

    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)

    PinCFlow_dev.calc_RHS(rhs, semi, 1.0)
    dt = 0.3
    errFlagBicg = false
    nIter = 0
    facprs = 1.0
    facray = 1.0
    PinCFlow_dev.poissonSolver(rhs, semi, dt, errFlagBicg, nIter, "expl", facray, facprs)

    test_arr(cache.dp, 6954.114011644017, 406.0728580408809, 32.1315910002854, tol = 1e-9)

    PinCFlow_dev.poissonSolver(rhs, semi, dt, errFlagBicg, nIter, "impl", facray, facprs)

    test_arr(cache.dp, 6949.851010678487, 405.8250898219831, 32.1214829208721, tol = 1e-9)

    PinCFlow_dev.pressureBoundaryCondition(semi)

    test_arr(cache.dp, 26794.367327608583, 799.5740338430224, 32.1214829208721, tol = 1e-9)

    PinCFlow_dev.correctorStep(semi, dt, "expl", facray, facprs)

    test_arr(var.exner, 26794.367327608583, 799.5740338430224, 32.1214829208721, tol = 1e-9)
    test_arr(var.u, 2996.1961128135836, 250.73352464698718, 64.47576033368075, tol = 1e-9)

    PinCFlow_dev.correctorStep(semi, dt, "impl", facray, facprs)

    set_lin_array_normalizer!.((kr_sp_tfc, kr_sp_w_tfc))
    @show get_norms(kr_sp_tfc)
    @show get_norms(kr_sp_w_tfc)
    test_arr(var.exner, 53588.734655217166, 1599.148067686045, 64.2429658417442, tol = 1e-9)
    test_arr(kr_sp_tfc, 368.64, 11.368166079012056, 0.54, tol = 1e-9)
    test_arr(kr_sp_w_tfc, 368.64, 11.368166079012056, 0.54, tol = 1e-9)
    test_arr(var.u, 1534.76, 23.24188460516906, 0.6, tol = 1e-9)
    test_arr(cache.corX, 1643.1243369208187, 248.88860795404267, 64.05576033368075,
             tol = 1e-9)
    test_arr(cache.corY, 0.15327530873362505, 0.013421375877617665, 0.0033422640397776215,
             tol = 1e-9)
    test_arr(var.exner, 53588.734655217166, 1599.148067686045, 64.2429658417442, tol = 1e-9)

    PinCFlow_dev.bicgstab(rhs, dt, semi, cache.sol_bicg, nIter, errFlagBicg, "expl")

    @show get_norms(cache.sol_bicg)
    @show get_norms(cache.p_bicg)
    test_arr(cache.sol_bicg, 1142.0651583477945, 66.0122418332421, 4.178712816389702,
             tol = 1e-9)
    test_arr(cache.p_bicg, 5.676571678854043e-6, 3.743420434190437e-7,
             3.6087341096433547e-8,
             tol = 1e-9)

    PinCFlow_dev.bicgstab(rhs, dt, semi, cache.sol_bicg, nIter, errFlagBicg, "impl")

    test_arr(cache.sol_bicg, 1142.0651583477945, 66.0122418332421, 4.178712816389702,
             tol = 1e-9)
    test_arr(cache.p_bicg, 5.676571678854043e-6, 3.743420434190437e-7,
             3.6087341096433547e-8,
             tol = 1e-9)
end
