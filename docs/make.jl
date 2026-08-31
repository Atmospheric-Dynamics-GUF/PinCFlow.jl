using Changelog: Changelog
using Documenter
using CairoMakie
using Revise
using PinCFlow

# Insert the example functions.
for folder in ("src/Examples/", "src/Examples/WavePacketTools/")
    for script_file in readdir(folder)
        if endswith(script_file, ".jl")
            script = read(folder * script_file, String)
            code = Regex(
                "(?sm)(?<=^```julia\\n)# " *
                folder *
                script_file *
                "(.(?!^```\\n))*",
            )
            if script_file == "periodic_hill.jl"
                page_files = ("README.md", "docs/src/examples.md")
            else
                page_files = ("docs/src/examples.md",)
            end
            for page_file in page_files
                if isfile(page_file)
                    page = replace(read(page_file, String), code => script)
                    open(page_file, "w") do io
                        write(io, page)
                        return
                    end
                end
            end
        end
    end
end

# Create the example plots.
mkpath("docs/src/examples/results/")
for name in names(PinCFlow.Examples)
    example = getproperty(PinCFlow.Examples, name)
    if example isa Function
        mktempdir() do directory
            example(;
                display_figure = false,
                output_file = directory * "/$(example).h5",
                plot_file = "docs/src/examples/results/$(example).svg",
            )
            return
        end
    end
end

# Copy the README file and use it as landing page of the docs.
cp(joinpath(dirname(@__DIR__), "README.md"), "docs/src/index.md"; force = true)

# Copy some files from the repository root directory to the docs and modify them
# as necessary.

open(joinpath(@__DIR__, "src", "authors.md"), "w") do io
    println(
        io,
        """
        ```@meta
        EditURL = "https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/blob/main/AUTHORS.md"
        ```
        """,
    )
    for line in eachline(joinpath(dirname(@__DIR__), "AUTHORS.md"))
        line = replace(
            line,
            "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)",
        )
        println(io, line)
    end
end

open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    println(
        io,
        """
        ```@meta
        EditURL = "https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/blob/main/CODE_OF_CONDUCT.md"
        ```
        """,
    )
    println(io, "# Code of Conduct")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, "  > ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    println(
        io,
        """
        ```@meta
        EditURL = "https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/blob/main/CONTRIBUTING.md"
        ```
        """,
    )
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, line)
    end
end

open(joinpath(@__DIR__, "src", "license.md"), "w") do io
    println(
        io,
        """
        ```@meta
        EditURL = "https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/blob/main/LICENSE.md"
        ```
        """,
    )
    println(io, "# License")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "LICENSE.md"))
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, "  > ", line)
    end
end

# Create a changelog.
Changelog.generate(
    Changelog.Documenter(),                        # output type
    joinpath(@__DIR__, "..", "NEWS.md"),           # input file
    joinpath(@__DIR__, "src", "changelog_tmp.md"); # output file
    repo = "Atmospheric-Dynamics-GUF/PinCFlow.jl",
    branch = "main",
)

# Fix the edit URL of the changelog.
open(joinpath(@__DIR__, "src", "changelog.md"), "w") do io
    for line in eachline(joinpath(@__DIR__, "src", "changelog_tmp.md"))
        if startswith(line, "EditURL")
            line = "EditURL = \"https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/blob/main/NEWS.md\""
        end
        println(io, line)
    end
end

# Remove the temporary file.
rm(joinpath(@__DIR__, "src", "changelog_tmp.md"))

# Generate the documentation.
makedocs(;
    sitename = "PinCFlow.jl",
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Theory" => [
            "Physics" => "theory/physics.md",
            "Numerics" => "theory/numerics.md",
        ],
        "Public API" => "public_api.md",
        "Internal reference" => [
            "PinCFlow" => "internal_reference/pincflow.md",
            "Macros" => "internal_reference/macros.md",
            "Types" => "internal_reference/types.md",
            "MPIOperations" => "internal_reference/mpi_operations.md",
            "Boundaries" => "internal_reference/boundaries.md",
            "FluxCalculator" => "internal_reference/flux_calculator.md",
            "PoissonSolver" => "internal_reference/poisson_solver.md",
            "Update" => "internal_reference/update.md",
            "MSGWaM" => "internal_reference/msgwam.md",
            "Integration" => "internal_reference/integration.md",
            "Output" => "internal_reference/output.md",
            "PinCFlowMakieExt" => "internal_reference/pincflow_makie_ext.md",
        ],
        "Developer guide" => "developer_guide.md",
        "Changelog" => "changelog.md",
        "Authors" => "authors.md",
        "Contributing" => "contributing.md",
        "Code of Conduct" => "code_of_conduct.md",
        "License" => "license.md",
    ],
    pagesonly = true,
    format = Documenter.HTML(;
        repolink = "https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl.git",
        size_threshold = nothing,
        size_threshold_warn = nothing,
    ),
)

# Only push previews if all the relevant environment variables are non-empty.
deploydocs(;
    repo = "github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl",
    devbranch = "main",
    push_preview = all(
        !isempty,
        (get(ENV, "GITHUB_TOKEN", ""), get(ENV, "DOCUMENTER_KEY", "")),
    ),
)
