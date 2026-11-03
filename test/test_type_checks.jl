function test_type_checks()
    @testset "Type checks" begin
        count = 0
        test_files() do root, file
            path = abspath(joinpath(root, file))
            content = Meta.parseall(read(path, String))
            postwalk(content) do x
                if @capture(x, (typeof(a_) <: b_) || (b_ >: typeof(a_)))
                    println("Bad type checks detected in the expression")
                    println(x)
                    println("of the file $(path)!")
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
