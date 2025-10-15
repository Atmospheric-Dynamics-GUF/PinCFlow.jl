using Pkg

Pkg.activate("format")

using JuliaFormatter
using Revise

format(".")

for (root, directories, files) in walkdir(".")
    for file in files
        if endswith(joinpath(root, file), r"\.(jl|md)")
            script = read(joinpath(root, file), String)

            for code_block in
                eachmatch(r"(?s)(?<=\n`{3}julia\n)(.(?!\n`{3}))*.", script)
                script = replace(
                    script,
                    code_block.match => format_text(
                        string(code_block.match);
                        margin = 80,
                        always_for_in = true,
                        whitespace_typedefs = true,
                        whitespace_ops_in_indices = true,
                        remove_extra_newlines = true,
                        pipe_to_function_call = true,
                        short_to_long_function_def = true,
                        force_long_function_def = true,
                        always_use_return = true,
                        normalize_line_endings = "unix",
                        separate_kwargs_with_semicolon = true,
                    ),
                )
            end

            open(joinpath(root, file), "w") do io
                write(io, script)
                return
            end
        end
    end
end
