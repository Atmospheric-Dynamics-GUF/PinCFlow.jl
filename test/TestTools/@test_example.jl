"""
```julia
@test_example(file::AbstractString, args::Vararg{Any})
```

Test the example defined in the run script `file`, using TrixiTest.

# Arguments

  - `file`: Example run script.

  - `args`: Keyword arguments that configure the test.
"""
macro test_example end

macro test_example(file::AbstractString, args::Vararg{Any})
    local l2 = get_kwarg(args, :l2, nothing)
    local linf = get_kwarg(args, :linf, nothing)
    local RealT_symbol = get_kwarg(args, :RealT, :Float64)
    RealT = getfield(@__MODULE__, RealT_symbol)
    atol_default = 500 * eps(RealT)
    rtol_default = sqrt(eps(RealT))
    local atol = get_kwarg(args, :atol, atol_default)
    local rtol = get_kwarg(args, :rtol, rtol_default)

    local kwargs = Pair{Symbol, Any}[]
    for arg in args
        if (
            arg.head == :(=) && !(
                arg.args[1] in
                (:additional_ignore_content, :l2, :linf, :RealT, :atol, :rtol)
            )
        )
            push!(kwargs, Pair(arg.args...))
        end
    end

    quote
        trixi_include(@__MODULE__, $(esc(file)); $kwargs...)

        l2_measured, linf_measured = invokelatest((@__MODULE__).compute_norms)
        @test length($l2) == length(l2_measured)
        for (l2_expected, l2_actual) in zip($l2, l2_measured)
            @test isapprox(l2_expected, l2_actual, atol = $atol, rtol = $rtol)
        end

        @test length($linf) == length(linf_measured)
        for (linf_expected, linf_actual) in zip($linf, linf_measured)
            @test isapprox(
                linf_expected,
                linf_actual,
                atol = $atol,
                rtol = $rtol,
            )
        end
    end
end
