using NCDatasets, PythonPlot

# adjust path to your netcdf file
ncpath = "/home/dolaptch/PF/runs/tjl04/pincflow_data_out.nc"
ds = Dataset(ncpath, "r")

# read vars (example: group "icevar", variable "sn")
sn = ds["icevar/sn"][:]        # dims typically (nx, ny, nz, time)
x  = ds["x"][:] .* 0.001       # km
z  = ds["z"][:, :, :] .* 0.001 # km (adjust if z is 1D)

tidx = size(sn, 4)             # last time index
fld = sn[:, 1, :, tidx]       # slice at iy=1 -> (nx, nz)

# If x and z are 1D, you can call pcolormesh(x, z, fld')
# If z is 2D use matching shapes; adapt below to your dimensions:
pcolormesh(x, reshape(z[:,1, :], size(z,1), size(z,3)), fld'; cmap="viridis")
colorbar()
title("sn (iy=1, t=$(tidx))")
savefig("examples/visualization/sn_pythonplot.png")

close(ds)