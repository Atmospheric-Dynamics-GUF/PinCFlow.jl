function test_ascii()
    @testset "ASCII" begin
        count = 0
        test_files() do root, file
            path = abspath(joinpath(root, file))
            content = read(path, String)
            non_ascii = r"[^\x00-\x7F]"
            if occursin(non_ascii, content)
                println("Non-ASCII characters")
                println(Tuple(m for m in eachmatch(non_ascii, content)))
                println("detected in $(path)!")
                println("")
                count += 1
            end
        end
        @test count == 0
    end

    nothing
end
