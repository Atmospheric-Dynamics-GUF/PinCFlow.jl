#!/usr/bin/env julia
# Format all Julia files in the project
# Usage: julia --project=format format.jl [path]
# If no path is provided, formats the entire project

using JuliaFormatter

path = length(ARGS) > 0 ? ARGS[1] : "."
println("Formatting: $path")
format(path)
