module Turbulence_research

using Revise
using PinCFlow
using MPI
using HDF5

include("wp_1d.jl")
include("wkb_wp_1d.jl")
include("wkb_wp.jl")
include("wp.jl")

export wp_1d, wkb_wp_1d, wkb_wp, wp

end
