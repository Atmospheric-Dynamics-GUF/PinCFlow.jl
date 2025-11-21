"""
```julia
test_example(
    file::AbstractString,
    reference::NTuple{2, <:NamedTuple},
    assignments::Vararg{Pair{Symbol, <:Any}};
    update_references::Bool = false,
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

  - `update_references`: Switch for updating the references in the test scripts instead of testing against them.

  - `atol`: Absolute tolerance that is passed as keyword argument to `isapprox`.

  - `rtol`: Relative tolerance that is passed as keyword argument to `isapprox`.
"""
function test_example end

function test_example(
    file::AbstractString,
    reference::NTuple{2, <:NamedTuple},
    assignments::Vararg{Pair{Symbol, <:Any}};
    update_references::Bool = false,
    atol::Real = 0,
    rtol::Real = 0,
)
    # Read the example script and modify it.
    script = read(file, String)
    start = findfirst(r"(?m)^ *integrate\(", script)[1]
    stop = Meta.parse(script, start)[2] - 1
    script = script[1:stop]
    script = replace(
        script,
        r"(?m)^ *using +Pkg *\n+" => "",
        r"(?m)^ *Pkg.activate\( *\"examples\" *\) *\n+" => "",
        r"(?m)^ *using +HDF5 *\n+" => "",
        r"(?m)^ *using +CairoMakie *\n+" => "",
    )
    script = replace_assignments(script, assignments...)

    # Run the modified script.
    eval(Meta.parseall(script))

    # Get the norms.
    (l2ref, linfref) = reference
    (l2, linf) = compute_norms()

    # Update the references or test against them.
    if update_references
        test_file = splitpath(file)[end]
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
