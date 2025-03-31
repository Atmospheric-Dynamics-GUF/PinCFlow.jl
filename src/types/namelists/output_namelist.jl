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
  A <: AbstractVector{<:AbstractVariable},
  B <: Bool,
  C <: Integer,
  D <: AbstractFloat,
}
  atmvarout::A
  prepare_restart::B
  restart::B
  iin::C
  output_steps::B
  noutput::C
  maxiter::C
  outputtimediff::D
  maxtime::D
  fancy_namelists::B
end

function OutputNamelist(;
  atmvarout = Vector{AbstractVariable}(undef, 0),
  prepare_restart = false,
  restart = false,
  iin = -1,
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
    output_steps,
    noutput,
    maxiter,
    outputtimediff,
    maxtime,
    fancy_namelists,
  )
end
