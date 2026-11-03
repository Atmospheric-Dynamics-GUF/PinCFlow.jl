function test_files(test_function::Function)
    for (root, directories, files) in walkdir(joinpath(@__DIR__, "../src"))
        for file in files
            test_function(root, file)
        end
    end

    nothing
end
