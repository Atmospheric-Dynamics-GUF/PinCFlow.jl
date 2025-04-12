module Output

using Dates
using HDF5
using ..Types
using ..Update

include("create_output.jl")
include("read_input!.jl")
include("write_output.jl")

export create_output, read_input!, write_output

end