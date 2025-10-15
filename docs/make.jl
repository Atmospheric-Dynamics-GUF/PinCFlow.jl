using Pkg

Pkg.activate("docs")

using Changelog: Changelog
using Documenter
using Revise
using PinCFlow

# Insert the example scripts.
@ivy for folder in ("examples/submit/", "examples/visualization/"),
    script_file in readdir(folder)

    if endswith(script_file, ".jl")
        script = read(folder * script_file, String)
        code = Regex(
            "(?s)(?<=\\n`{3}julia\\n)# " *
            folder *
            script_file[1:(end - 3)] *
            "\\.jl\\n(.(?!\\n`{3}))*.",
        )
        if script_file == "periodic_hill.jl"
            page_file = "README.md"
        else
            page_file =
                "docs/src/examples/" *
                script_file[1:(end - 3)] *
                "_simulation.md"
        end
        page = replace(read(page_file, String), code => script)
        open(page_file, "w") do io
            write(io, page)
            return
        end
    end
end

# Copy the example plots.
mkpath("docs/src/examples/results/")
for file in readdir("examples/results/"; join = true)
    cp(file, "docs/src/" * file; force = true)
end

# Copy the README file and use it as landing page of the docs.
cp(joinpath(dirname(@__DIR__), "README.md"), "docs/src/index.md"; force = true)

# Copy list of authors to not need to synchronize it manually.
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(
    authors_text,
    "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)",
)
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Copy some files from the repository root directory to the docs and modify them
# as necessary.

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
        "Examples" => [
            "Mountain-wave simulation" => "examples/mountain_wave_simulation.md",
            "WKB mountain-wave simulation" => "examples/wkb_mountain_wave_simulation.md",
        ],
        "Theory" => [
            "Physics" => "theory/physics.md",
            "Numerics" => "theory/numerics.md",
        ],
        "Reference" => [
            "PinCFlow" => "reference/pincflow.md",
            "Types" => "reference/types.md",
            "MPIOperations" => "reference/mpi_operations.md",
            "Boundaries" => "reference/boundaries.md",
            "FluxCalculator" => "reference/flux_calculator.md",
            "PoissonSolver" => "reference/poisson_solver.md",
            "Update" => "reference/update.md",
            "MSGWaM" => "reference/msgwam.md",
            "Integration" => "reference/integration.md",
            "Output" => "reference/output.md",
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

deploydocs(;
    repo = "github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl",
    devbranch = "main",
    push_preview = true,
)
