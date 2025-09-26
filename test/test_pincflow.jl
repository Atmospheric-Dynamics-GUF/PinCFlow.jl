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
