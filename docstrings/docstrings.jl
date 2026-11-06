using CairoMakie
using PinCFlow

using .Meta: parse, parseall
using MacroTools: postwalk
using REPL: moduleusings

PinCFlowMakieExt = Base.get_extension(PinCFlow, :PinCFlowMakieExt)

for (folder, module_name) in (("src", "PinCFlow"), ("ext", "PinCFlowMakieExt"))
    for (root, directories, files) in walkdir(folder)
        for file in files
            path = abspath(joinpath(root, file))
            content_string = read(path, String)
            content_expression = parseall(content_string)
            module_file = false
            postwalk(content_expression) do x
                x isa Expr && x.head === :module && (module_file = true)
                x
            end

            parent_name = replace(
                module_name * rsplit(dirname(path), "/$folder")[end],
                "/" => ".",
            )
            object_name = basename(path)[1:(end - 3)]

            parent = eval(Meta.parse(parent_name))
            object = getproperty(parent, Symbol(object_name))

            if module_file
                scopes = sort([
                    "$scope" for scope in moduleusings(object) if
                    startswith("$scope", "PinCFlow")
                ])
                exports = sort([
                    "$(parentmodule(getproperty(object, name))).$name" for
                    name in names(object) if "$name" != object_name
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
                            call = Core.eval(parent, x.args[1])
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

            !(x isa String) && continue

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

                    docstring =
                        "\"\"\"\n" *
                        strip(
                            occursin(pattern, docstring) ?
                            replace(docstring, pattern => section) :
                            rstrip(docstring, '\n') * section,
                            '\n',
                        ) *
                        "\n\"\"\""

                    content_string =
                        content_string[1:(start - 1)] *
                        prefix *
                        docstring *
                        content_string[(stop + 1):end]

                    stop = start - 1 + length(prefix) + length(docstring)

                    open(path, "w") do io
                        write(io, content_string)
                    end
                end
            end
        end
    end
end
