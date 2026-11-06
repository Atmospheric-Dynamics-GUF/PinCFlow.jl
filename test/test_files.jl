function test_files(test_function::Function)
    for folder in (joinpath(@__DIR__, "../src"), joinpath(@__DIR__, "../ext"))
        for (root, directories, files) in walkdir(folder)
            for file in files
                test_function(root, file)
            end
        end
    end

    nothing
end
