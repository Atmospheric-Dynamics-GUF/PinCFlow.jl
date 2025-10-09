# Testing 
using Test
using LinearAlgebra: norm
using TrixiTest: trixi_include, get_kwarg, append_to_kwargs
using HDF5
import PinCFlow: PinCFlow

function examples_dir()
    return pkgdir(PinCFlow, "examples/submit")
end

submit_dir = examples_dir()

macro test_example(file, args...)
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
        local file_to_include = $(esc(file))
        local content = read(file_to_include, String)
        content = replace(content, r"Pkg\.activate\(.*?\)\n?" => "")
        local temp_file = tempname() * ".jl"
        write(temp_file, content)

        trixi_include(@__MODULE__, temp_file; $kwargs...)
        rm(temp_file; force = true)

        l2_measured, linf_measured =
            invokelatest((@__MODULE__).compute_norms, "./pincflow_output.h5")
        #        l2_measured, linf_measured = compute_norms("./pincflow_output.h5")        
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
#TODO: filepath is not used
function compute_norms(filepath::String)
    h5open("./pincflow_output.h5", "r") do data
        vars = [k for k in keys(data) if !(k in ["x", "y", "z"])]
        l2 = [norm(read(data[v]), 2) for v in vars]
        linf = [norm(read(data[v]), Inf) for v in vars]
        return l2, linf
    end
end
