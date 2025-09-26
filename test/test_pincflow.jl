# Testing 
using Test
using LinearAlgebra: norm

# Used for testing
function examples_dir()
    return pkgdir(PinCFlow, "examples/submit")
end

export examples_dir

function test_example(
    file::AbstractString,
    assignments::AbstractString...;
    l2 = nothing,
    linf = nothing,
    atol = 500 * eps(Float32),
    rtol = sqrt(eps(Float32)),
)
    script = read(file, String)

    for assignment in assignments
        var_name, new_value = split(assignment, "="; limit = 2)
        var_name = strip(var_name)
        new_value = strip(new_value)

        pattern = Regex("\\b" * var_name * "\\s*=\\s*[^,\\)\\n]+")
        replacement = var_name * " = " * new_value
        script = replace(script, pattern => replacement)
    end

    eval(Meta.parseall(script))

    l2_measured, linf_measured = compute_norms("./pincflow_output.h5")

    @test length(l2) == length(l2_measured)
    for (l2_expected, l2_actual) in zip(l2, l2_measured)
        @test isapprox(l2_expected, l2_actual, atol = atol, rtol = rtol)
    end

    @test length(linf) == length(linf_measured)
    for (linf_expected, linf_actual) in zip(linf, linf_measured)
        @test isapprox(linf_expected, linf_actual, atol = atol, rtol = rtol)
    end
end

function compute_norms(filepath::String)
    h5open("./pincflow_output.h5", "r") do data
        vars = [k for k in keys(data) if k ∉ ["x", "y", "z"]]
        l2 = [norm(read(data[v]), 2) for v in vars]
        linf = [norm(read(data[v]), Inf) for v in vars]
        return l2, linf
    end
end
