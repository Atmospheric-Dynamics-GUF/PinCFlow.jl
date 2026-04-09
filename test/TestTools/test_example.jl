"""
```julia
test_example(
    example::Function,
    keywords::NamedTuple,
    reference::NTuple{2, <:NamedTuple};
    update_references::Bool = false,
    atol::Real = 0,
    rtol::Real = 0,
)
```

Run an example simulation with `keywords`, compute the ``L_2`` and ``L_\\infty`` norms of each output variable and test them against the `reference`.

# Arguments

  - `example`: Example simulation function.

  - `keywords`: Keyword arguments to pass to `example`.

  - `reference`: Reference ``L_2`` and ``L_\\infty`` norms.

# Keywords

  - `update_references`: Switch for updating the references in the test scripts instead of testing against them.

  - `atol`: Absolute tolerance that is passed as keyword argument to `isapprox`.

  - `rtol`: Relative tolerance that is passed as keyword argument to `isapprox`.
"""
function test_example end

function test_example(
    example::Function,
    keywords::NamedTuple,
    reference::NTuple{2, <:NamedTuple};
    update_references::Bool = false,
    atol::Real = 0,
    rtol::Real = 0,
)

    # Call the example function with the provided keywords.
    redirect_stdio(; stderr = devnull, stdout = devnull) do
        example(; keywords...)
        return
    end

    # Get the norms.
    (l2ref, linfref) = reference
    (l2, linf) = compute_norms()

    # Update the references or test against them.
    if update_references
        test_file = "test_" * string(nameof(example)) * ".jl"
        script = replace_assignments(
            read(test_file, String),
            :l2 => l2,
            :linf => linf,
        )
        open(test_file, "w") do io
            write(io, script)
            return
        end
        format(test_file)
    else
        for (l, lref, label) in
            ((l2, l2ref, "L2 norms"), (linf, linfref, "LInf norms"))
            @testset "$label" begin
                for key in keys(lref)
                    @testset "$key" begin
                        @test isapprox(
                            l[key],
                            lref[key];
                            atol,
                            rtol = rtol == 0 ?
                                   (atol > 0 ? 0 : sqrt(eps(Float32))) : rtol,
                        )
                    end
                end
            end
        end
    end

    return
end
