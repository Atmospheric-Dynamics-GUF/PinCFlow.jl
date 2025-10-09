"""
```julia
test_example(
    file::AbstractString,
    reference::NTuple{2, <:NamedTuple},
    assignments::Vararg{Pair{Symbol, <:Any}};
    test::Bool = true,
    atol::Real = 0,
    rtol::Real = 0,
)
```

Run an example simulation with modified `assignments`, compute the ``L_2`` and ``L_\\infty`` norms of each output variable and test them against the `reference`.

# Arguments

  - `file`: Example run script.

  - `reference`: Reference ``L_2`` and ``L_\\infty`` norms.

  - `assignments`: Replacements for assignments in the example run script.

# Keywords

  - `test`: Test the ``L_2`` and ``L_\\infty`` norms against the reference if set to `true`, print them otherwise.

  - `atol`: Absolute tolerance that is passed as keyword argument to `isapprox`.

  - `rtol`: Relative tolerance that is passed as keyword argument to `isapprox`.
"""
function test_example end

function test_example(
    file::AbstractString,
    reference::NTuple{2, <:NamedTuple},
    assignments::Vararg{Pair{Symbol, <:Any}};
    test::Bool = true,
    atol::Real = 0,
    rtol::Real = 0,
)
    # Read the example script and replace assignments.
    script = replace(
        read(file, String),
        r"(?m)^ *using +Pkg *\n+" => "",
        r"(?m)^ *Pkg.activate\( *\"examples\" *\) *\n+" => "",
    )
    for assignment in assignments
        (name, value) = assignment
        range = findfirst(Regex("$name *= *"), script)
        if range !== nothing
            (start, stop) = extrema(range)
            suffix = ""
            stop += 1
            try
                stop = Meta.parse(script, stop)[2]
                suffix = "\n"
            catch
                while !(script[stop] in (',', ')'))
                    stop = Meta.parse(script, stop; greedy = false)[2]
                end
            end
            stop -= 1
            script =
                replace(script, script[start:stop] => "$name = $value" * suffix)
        end
    end

    # Run the modified script.
    eval(Meta.parseall(script))

    # Get the norms.
    (l2ref, linfref) = reference
    (l2, linf) = compute_norms()

    # Test or print the norms.
    if test
        for (l, lref, label) in
            ((l2, l2ref, "L2 norms"), (linf, linfref, "LInf norms"))
            @testset "$label" begin
                for key in keys(lref)
                    @testset "$key" begin
                        @test isapprox(l[key], lref[key]; atol, rtol)
                    end
                end
            end
        end
    else
        println("l2 = ", l2)
        println("linf = ", linf)
    end

    return
end
