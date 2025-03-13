abstract type AbstractOutput end
struct RHOP <: AbstractOutput end
struct U <: AbstractOutput end
struct US <: AbstractOutput end
struct V <: AbstractOutput end
struct VS <: AbstractOutput end
struct W <: AbstractOutput end
struct WS <: AbstractOutput end
struct WTFC <: AbstractOutput end
struct WSTFC <: AbstractOutput end
struct THETAP <: AbstractOutput end
struct PIP <: AbstractOutput end

struct OutputNamelist{
    A<:Vector{<:AbstractOutput},
    B<:Bool,
    C<:Integer,
    D<:String,
    E<:AbstractFloat,
}
    atmvarout::A
    prepare_restart::B
    restart::B
    iin::C
    runname::D
    output_steps::B
    noutput::C
    maxiter::C
    outputtimediff::E
    maxtime::E
    fancy_namelists::B
end

function OutputNamelist(;
  atmvarout = Vector{AbstractOutput}(undef, 0),
  prepare_restart = false,
  restart = false,
  iin = -1,
  runname = "",
  output_steps = false,
  noutput = 1,
  maxiter = 1,
  outputtimediff = 3600.0,
  maxtime = 3600.0,
  fancy_namelists = true,
)
  return OutputNamelist(
    atmvarout,
    prepare_restart,
    restart,
    iin,
    runname,
    output_steps,
    noutput,
    maxiter,
    outputtimediff,
    maxtime,
    fancy_namelists,
  )
end
