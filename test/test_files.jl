function test_files(
    test_function::Function;
    folder::AbstractString = abspath(joinpath(@__DIR__, "../src")),
)
    for (root, directories, files) in walkdir(folder)
        for file in files
            test_function(root, file)
        end
    end

    nothing
end
