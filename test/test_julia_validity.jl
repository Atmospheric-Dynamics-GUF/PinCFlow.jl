function test_julia_validity()
    @testset "Julia validity" begin
        count = 0
        test_files() do root, file
            path = abspath(joinpath(root, file))
            (name, extension) = splitext(file)
            content = Meta.parseall(read(path, String))
            if !(
                extension == ".jl" &&
                content isa Expr &&
                content.head === :toplevel
            )
                println("$path is not a valid Julia file!")
                count += 1
            end
        end
        @test count == 0
    end

    nothing
end
