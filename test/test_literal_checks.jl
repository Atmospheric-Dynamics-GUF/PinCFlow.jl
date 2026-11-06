function test_literal_checks()
    @testset "Literal checks" begin
        count = 0
        test_files() do root, file
            path = abspath(joinpath(root, file))
            content = Meta.parseall(read(path, String))
            postwalk(content) do x
                if (@capture(x, a_ == b_) || @capture(x, a_ != b_)) && (
                    (@capture(a, Val(c_quote)) && c.args[1] isa Symbol) ||
                    (@capture(b, Val(c_quote)) && c.args[1] isa Symbol) ||
                    (a isa QuoteNode && a.value isa Symbol) ||
                    (b isa QuoteNode && b.value isa Symbol) ||
                    a === :(Val(nothing)) ||
                    b === :(Val(nothing)) ||
                    a === :nothing ||
                    b === :nothing
                )
                    println("Bad literal check")
                    println(x)
                    println("detected in $(path)!")
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
