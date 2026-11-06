function test_docstrings(;
    fix::Bool = false,
    folder::AbstractString = abspath(joinpath(@__DIR__, "../src")),
)
    if abspath(folder) == abspath(joinpath(@__DIR__, "../src"))
        module_name = "PinCFlow"
    elseif abspath(folder) == abspath(joinpath(@__DIR__, "../ext"))
        module_name = "PinCFlowMakieExt"
    else
        error("Unknown folder!")
    end

    folder = abspath(folder)

    up_to_date = true

    test_files(; folder) do root, file
        path = abspath(joinpath(root, file))
        content_string = read(path, String)
        content_expression = parseall(content_string)
        module_file = false
        postwalk(content_expression) do x
            x isa Expr && x.head === :module && (module_file = true)
            x
        end

        location = eval(
            parse(
                replace(
                    module_name * split(dirname(path), "$folder")[end],
                    "/" => ".",
                ),
            ),
        )
        object = getproperty(location, Symbol(basename(path)[1:(end - 3)]))
        parent = parentmodule(object)

        if module_file
            scopes = sort([
                "$scope" for scope in moduleusings(object) if
                startswith("$scope", "PinCFlow")
            ])
            exports = sort([
                "$(parentmodule(getproperty(object, name))).$name" for
                name in names(object) if name !== nameof(object)
            ])
            list_names = ["Scopes", "Exports"]
            lists = [scopes, exports]
        else
            calls = Set{Symbol}()
            postwalk(content_expression) do x
                x isa Expr &&
                    x.head === :call &&
                    x.args[1] isa Symbol &&
                    try
                        call = Core.eval(location, x.args[1])
                        parent_module = parentmodule(call)
                        startswith("$parent_module", "PinCFlow") &&
                            call !== object &&
                            push!(
                                calls,
                                Symbol(parent_module, :., nameof(call)),
                            )
                    catch
                    end
                x
            end
            list_names = ["Calls"]
            lists = [sort(collect(calls))]
        end

        (x, start, stop) = (nothing, 1, 1)
        while stop < length(content_string)
            start = stop
            (x, stop) = parse(content_string, start)

            x isa Expr &&
                x.head === :macrocall &&
                x.args[1] isa GlobalRef &&
                x.args[1].mod === Core &&
                x.args[1].name === Symbol("@doc") &&
                break
        end
        (x, stop) = parse(content_string, start; greedy = false)
        stop -= 1

        if x isa String
            for (list_name, list) in zip(list_names, lists)
                (prefix, docstring) =
                    split(content_string[start:stop], r"\"+"; limit = 2)
                docstring = rstrip(docstring, '"')

                if length(list) > 0
                    pattern = Regex("\n+# $list_name(\n+[^#\n].+)*\n+")

                    section =
                        "\n\n# $list_name\n\n" *
                        join(("  - [`$item`](@ref)" for item in list), "\n\n") *
                        "\n\n"

                    new_docstring =
                        "\n" *
                        strip(
                            occursin(pattern, docstring) ?
                            replace(docstring, pattern => section) :
                            rstrip(docstring, '\n') * section,
                            '\n',
                        ) *
                        "\n"

                    if new_docstring != docstring
                        up_to_date = false

                        if fix
                            print(
                                "Updating the docstring of $parent.$(nameof(object))...",
                            )

                            bounded_docstring =
                                "\"\"\"" * new_docstring * "\"\"\""

                            content_string =
                                content_string[1:(start - 1)] *
                                prefix *
                                bounded_docstring *
                                content_string[(stop + 1):end]

                            stop =
                                start - 1 +
                                length(prefix) +
                                length(bounded_docstring)

                            open(path, "w") do io
                                write(io, content_string)
                            end

                            println("done.")
                        else
                            println(
                                "The docstring of $parent.$(nameof(object)) is not up to date!",
                            )
                        end
                    end
                end
            end
        end
    end

    !fix && @testset "Docstrings" begin
        @test up_to_date
    end
end
