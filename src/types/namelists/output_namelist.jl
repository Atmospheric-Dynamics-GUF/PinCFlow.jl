abstract type AbstractVariable end
struct Rho <: AbstractVariable end
struct RhoP <: AbstractVariable end
struct U <: AbstractVariable end
struct US <: AbstractVariable end
struct V <: AbstractVariable end
struct VS <: AbstractVariable end
struct W <: AbstractVariable end
struct WS <: AbstractVariable end
struct WTFC <: AbstractVariable end
struct WSTFC <: AbstractVariable end
struct ThetaP <: AbstractVariable end
struct PiP <: AbstractVariable end

struct OutputNamelist{
  A <: Vector{<:AbstractVariable},
  B <: Bool,
  C <: Integer,
  D <: String,
  E <: AbstractFloat,
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
  atmvarout = Vector{AbstractVariable}(undef, 0),
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
