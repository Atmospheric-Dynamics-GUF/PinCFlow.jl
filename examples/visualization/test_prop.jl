using Pkg

Pkg.activate("examples")

using HDF5
using CairoMakie
using Revise
using PinCFlow

h5open("wkb_wave_packet.h5", "r") do data; 
    x = data["x"][:] ./ 1000
    y = data["y"][:] ./ 1000
    z = data["z"][:, :, :] ./ 1000
    t = data["t"][:] ./ 3600
    kp = data["kp"][:]  .* 1000
    m = data["m"][:] .*1000
    (nt,) = size(t)
    (nx, ny, nz) = size(z)
    #x = [xi for xi in x, j in 1:ny, k in 1:nz, l in 1:nt]
    #y = [yj for i in 1:nx, yj in y, k in 1:nz]
    # Get the time.
    

    # Get the variable.
    phi = data["wavespectrum"][8, 8, 16,:, :, 3]
    #print(size(x[:, 1, 2, 2]))
    #print(size(y))
    #print(size(z))
    println(size(phi))
    println(phi) 


end


