function test_formatting(; directory::AbstractString = "..", fix::Bool = false)
    keywords = (
        always_for_in = true,
        always_use_return = true,
        force_long_function_def = true,
        margin = 80,
        normalize_line_endings = "unix",
        pipe_to_function_call = true,
        remove_extra_newlines = true,
        separate_kwargs_with_semicolon = true,
        short_to_long_function_def = true,
        whitespace_ops_in_indices = true,
        whitespace_typedefs = true,
    )

    formatted = true
    for (root, directories, files) in walkdir(directory)
        for file in files
            if endswith(joinpath(root, file), r"\.(jl|md)")
                script = read(joinpath(root, file), String)
                modified_script = script

                # Remove trailing whitespaces.
                modified_script = join(
                    (rstrip(line) for line in split(modified_script, "\n")),
                    "\n",
                )

                # Format Julia-code blocks in docstrings.
                for code_block in
                    eachmatch(r"(?sm)(?<=^```julia\n)(.(?!^```\n))*", script)
                    modified_script = replace(
                        modified_script,
                        code_block.match => format_text(
                            string(code_block.match);
                            keywords...,
                        ),
                    )
                end

                # Format everything else.
                if endswith(joinpath(root, file), ".jl")
                    modified_script = format_text(modified_script; keywords...)
                end

                # Update the file.
                if modified_script != script
                    formatted = false
                    fix && open(joinpath(root, file), "w") do io
                        formatted = false
                        write(io, modified_script)
                        return
                    end
                end
            end
        end
    end

    !fix && @testset "Formatting" begin
        @test formatted
    end

    return
end
