function test_long_named_function_definitions()
    @testset "Long named-function definitions" begin
        count = 0
        test_files() do root, file
            path = abspath(joinpath(root, file))
            content = Meta.parseall(read(path, String))
            postwalk(content) do x
                if isdef(x) &&
                   length(x.args) > 1 &&
                   :name in keys(splitdef(x)) &&
                   x.head !== :function
                    println("Short named-function definition")
                    println(x)
                    println("detected in $(path)!")
                    println("")
                    count += 1
                end
                x
            end
            nothing
        end
        @test count == 0
    end

    nothing
end
