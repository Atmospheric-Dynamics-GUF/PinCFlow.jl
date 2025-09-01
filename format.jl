using JuliaFormatter

for target in ("docs/make.jl", "examples/", "src/", "test/")
    format(target)
end
