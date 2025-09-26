# Testing 

using Test
using PinCFlow: examples_dir, analysis_callback
using TrixiTest: @test_trixi_include_base, append_to_kwargs

EXAMPLES_DIR = examples_dir()

macro test_pinc_include(expr, args...)
    local add_to_additional_ignore_content = [
        r"WARNING: Method definition .* overwritten .*.\n",
        r"Precompiling PinCFlow\.\.\.\n",
        r"   \d+\.\d+ ms  ✓ PinCFlow\n",
        r"  \d+ dependenc(y|ies) successfully precompiled in \d+ seconds\. \d+ already precompiled\.\n",
        r"Precompiling .+\.\.\.\n",
        r"   [\d\.]+ ms  ✓ .+\n",
        r"  \d+ dependenc(y|ies) successfully precompiled in [\d\.]+ seconds\. \d+ already precompiled\.\n",
    ]

    args = append_to_kwargs(
        args,
        :additional_ignore_content,
        add_to_additional_ignore_content,
    )

    quote
        @test_trixi_include_base($(esc(expr)), $(args...))
    end
end


function test_example(file::AbstractString, args...)

    local l2 = get_kwarg(args, :l2, nothing)
    local linf = get_kwarg(args, :linf, nothing)
    atol_default = 500 * eps(RealT)
    rtol_default = sqrt(eps(RealT))
    local atol = get_kwarg(args, :atol, atol_default)
    local rtol = get_kwarg(args, :rtol, rtol_default)

    script = read(file, String)
    modified_script = replace(
        script,
        r"\bsizex *= *\d+\b" => "sizex = 3",
        r"\bsizey *= *\d+\b" => "sizey = 3",
        r"\bsizez *= *\d+\b" => "sizez = 3",
        r"\bnpx *= *\d+\b" => "npx = 1",
        r"\bnpy *= *\d+\b" => "npy = 1",
        r"\bnpz *= *\d+\b" => "npz = 1",
    )
    eval(Meta.parseall(modified_script))

    l2_measured, linf_measured = invokelatest((@__MODULE__).PinCFlow.analysis_callback,
                                                      (@__MODULE__).sol)
	 @test length($l2) == length(l2_measured)
                for (l2_expected, l2_actual) in zip($l2, l2_measured)
                    @test isapprox(l2_expected, l2_actual, atol = $atol, rtol = $rtol)
                end
	
	@test length($linf) == length(linf_measured)
                for (linf_expected, linf_actual) in zip($linf, linf_measured)
                    @test isapprox(linf_expected, linf_actual, atol = $atol, rtol = $rtol)
                end

end

