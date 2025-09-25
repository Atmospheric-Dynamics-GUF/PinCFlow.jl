using Documenter
using PinCFlow

# Insert example scripts.
@ivy for folder in ("examples/submit/", "examples/visualization/"),
    script_file in readdir(folder)

    if endswith(script_file, ".jl") && script_file != "style.jl"
        script = read(folder * script_file, String)
        code = Regex(
            "(?s)(?<=`{3}julia\\n)# " *
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

# Copy example plots.
for file in readdir("examples/results/"; join = true)
    cp(file, "docs/src/" * file; force = true)
end

# Copy README file and SVGs.
cp("README.md", "docs/src/index.md"; force = true)
cp("pincflow_modules.svg", "docs/src/pincflow_modules.svg"; force = true)
cp("pincflow_structures.svg", "docs/src/pincflow_structures.svg"; force = true)

# Generate documentation.
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
    ],
    pagesonly = true,
    format = Documenter.HTML(;
        repolink = "https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl.git",
        size_threshold = nothing,
        size_threshold_warn = nothing,
    ),
)

deploydocs(repo = "github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl",
           devbranch = "main",
           push_preview = true)
