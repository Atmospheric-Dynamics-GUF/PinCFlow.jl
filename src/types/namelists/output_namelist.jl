struct OutputNamelist{
    A <: Tuple{Vararg{Symbol, <:Integer}},
    B <: Bool,
    C <: Integer,
    D <: AbstractFloat,
    E <: AbstractString,
}
    output_variables::A
    prepare_restart::B
    restart::B
    iin::C
    output_steps::B
    noutput::C
    maxiter::C
    outputtimediff::D
    maxtime::D
    fancy_namelists::B
    folder::E
end

function OutputNamelist(;
    output_variables = (),
    prepare_restart = false,
    restart = false,
    iin = -1,
    output_steps = false,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 3.6E+3,
    maxtime = 3.6E+3,
    fancy_namelists = true,
    folder = "./",
)
    return OutputNamelist(
        output_variables,
        prepare_restart,
        restart,
        iin,
        output_steps,
        noutput,
        maxiter,
        outputtimediff,
        maxtime,
        fancy_namelists,
        folder,
    )
end
