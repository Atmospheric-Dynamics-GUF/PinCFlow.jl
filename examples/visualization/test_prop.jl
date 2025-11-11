using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow
"""
h5open("prop_perturb.h5", "r") do data; 
    x = data["x"][:] ./ 1000
    y = data["y"][:] ./ 1000
    z = data["z"][:, :, :] ./ 1000
    t = data["t"][:] ./ 3600
    (nt,) = size(t)
    (nx, ny, nz) = size(z)
    #x = [xi for xi in x, j in 1:ny, k in 1:nz, l in 1:nt]
    #y = [yj for i in 1:nx, yj in y, k in 1:nz]
    # Get the time.
    

    # Get the variable.
    #phi = data["w"][1, 1, :, :]
    #print(size(x[:, 1, 2, 2]))
    #print(size(y))
    #print(size(z))
    #print(size(phi))
    #print(x) 


end
a = [0.1, 0.2]
b = [0.11 0.12]
c = [0.111, 0.112, 0.113]
a_array = [aj for aj in a, i in 1:2, j in 1:3]
b_array = [bj for i in 1:2, bj in b, j in 1:3]
c_array = [cj for i in 1:2, j in 1:2, cj in c]
print(b_array)
print((b_array[:, :, 1]))
"""
