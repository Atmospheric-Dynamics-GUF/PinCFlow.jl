"""
```julia
create_output(state::State)
```

Create an HDF5 output file with one dataset for each variable.

The dimensions of the datasets are set to those of the domain, whereas the chunks are set to the dimensions of the MPI subdomains, in preparation for parallel output. Datasets for the grid, i.e. the fields `x`, `y` and `zc` of `state.grid`, the time and the fields of `state.atmosphere` are always created, regardless of the specifications in `state.namelists.output`. The one exception to this is the Boussinesq mode, in which no datasets are created for the fields of `state.atmosphere`, since they do not have a spatial dependence.

# Arguments

  - `state`: Model state.
"""
function create_output end

function create_output(state::State)
	(; x_size, y_size, z_size, npx, npy, npz) = state.namelists.domain
	(; prepare_restart, save_ray_volumes, output_variables, output_file) =
		state.namelists.output
	(; model, test_case) = state.namelists.setting
	(; comm) = state.domain
	(; nray_max) = state.wkb
	(; nxnscxx, nynscyy, nznsczz) = state.ice.subgrid

	# Set the chunk dimensions.
	cr = nray_max
	cx = div(x_size, npx)
	cy = div(y_size, npy)
	cz = div(z_size, npz)
	ct = 1

	# Prepare the output.
	h5open(output_file, "w", comm) do file

		# Create datasets for the dimensions.
		create_dataset(file, "x", datatype(Float32), dataspace((x_size,)))
		create_dataset(file, "y", datatype(Float32), dataspace((y_size,)))
		create_dataset(
			file,
			"z",
			datatype(Float32),
			dataspace((x_size, y_size, z_size));
			chunk = (cx, cy, cz),
		)
		#changes for cloudcover
		if !(typeof(state.namelists.ice.icesetup) <: NoIce) && state.namelists.ice.cloudcover isa CloudCoverOn

			(; sizex2, sizey2, sizez2) = state.ice.subgrid

			cx2 = div(sizex2, npx)
			cy2 = div(sizey2, npy)
			cz2 = div(sizez2, npz)

			create_dataset(file, "x2", datatype(Float32), dataspace((sizex2,)))
			create_dataset(file, "y2", datatype(Float32), dataspace((sizey2,)))

			create_dataset(
				file,
				"z2",
				datatype(Float32),
				dataspace((sizex2, sizey2, sizez2));
				chunk = (cx2, cy2, cz2),
			)
		end

		create_dataset(
			file,
			"t",
			datatype(Float32),
			dataspace((0,), (-1,));
			chunk = (ct,),
		)

		# Create datasets for the background.
		if model != Boussinesq()
			for label in ("rhobar", "thetabar", "n2")
				create_dataset(
					file,
					label,
					datatype(Float32),
					dataspace((x_size, y_size, z_size));
					chunk = (cx, cy, cz),
				)
			end

			if model == Compressible()
				create_dataset(
					file,
					"p",
					datatype(Float32),
					dataspace(
						(x_size, y_size, z_size, 0),
						(x_size, y_size, z_size, -1),
					);
					chunk = (cx, cy, cz, ct),
				)
			else
				create_dataset(
					file,
					"p",
					datatype(Float32),
					dataspace((x_size, y_size, z_size));
					chunk = (cx, cy, cz),
				)
			end
		end

		# Create datasets for the prognostic variables.
		if prepare_restart || :rhop in output_variables
			create_dataset(
				file,
				"rhop",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if :u in output_variables
			create_dataset(
				file,
				"u",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if prepare_restart || :us in output_variables
			create_dataset(
				file,
				"us",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if :v in output_variables
			create_dataset(
				file,
				"v",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if prepare_restart || :vs in output_variables
			create_dataset(
				file,
				"vs",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if :w in output_variables
			create_dataset(
				file,
				"w",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if :ws in output_variables
			create_dataset(
				file,
				"ws",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if :wt in output_variables
			create_dataset(
				file,
				"wt",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if prepare_restart || :wts in output_variables
			create_dataset(
				file,
				"wts",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if :thetap in output_variables
			create_dataset(
				file,
				"thetap",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end
		if prepare_restart || :pip in output_variables
			create_dataset(
				file,
				"pip",
				datatype(Float32),
				dataspace(
					(x_size, y_size, z_size, 0),
					(x_size, y_size, z_size, -1),
				);
				chunk = (cx, cy, cz, ct),
			)
		end

		if !(typeof(state.namelists.tracer.tracer_setup) <: NoTracer)
			for field in fieldnames(TracerPredictands)
				create_dataset(
					file,
					string(field),
					datatype(Float32),
					dataspace(
						(x_size, y_size, z_size, 0),
						(x_size, y_size, z_size, -1),
					);
					chunk = (cx, cy, cz, ct),
				)
			end
		end

		if !(typeof(state.namelists.ice.icesetup) <: NoIce)
			for field in fieldnames(IcePredictands)
				create_dataset(
					file,
					string(field),
					datatype(Float32),
					dataspace(
						(x_size, y_size, z_size, 0),
						(x_size, y_size, z_size, -1),
					);
					chunk = (cx, cy, cz, ct),
				)
			end
		end

		if !(typeof(state.namelists.ice.icesetup) <: NoIce)
			for field in fieldnames(IceAuxiliaries)
				create_dataset(
					file,
					string(field),
					datatype(Float32),
					dataspace(
						(x_size, y_size, z_size, 0),
						(x_size, y_size, z_size, -1),
					);
					chunk = (cx, cy, cz, ct),
				)
			end
		end

		if state.namelists.tracer.leading_order_impact &&
		   :dchidt in output_variables
			for field in fieldnames(TracerWKBImpact)
				create_dataset(
					file,
					string(field),
					datatype(Float32),
					dataspace(
						(x_size, y_size, z_size, 0),
						(x_size, y_size, z_size, -1),
					);
					chunk = (cx, cy, cz, ct),
				)
			end
		end

		# Create datasets for WKB variables.
		if typeof(test_case) <: AbstractWKBTestCase

			# Create datasets for ray-volume properties.
			if prepare_restart || save_ray_volumes
				for field in (
					"xr",
					"yr",
					"zr",
					"dxr",
					"dyr",
					"dzr",
					"kr",
					"lr",
					"mr",
					"dkr",
					"dlr",
					"dmr",
					"nr",
				)
					create_dataset(
						file,
						field,
						datatype(Float32),
						dataspace(
							(nray_max, x_size, y_size, z_size + 1, 0),
							(nray_max, x_size, y_size, z_size + 1, -1),
						);
						chunk = (cr, cx, cy, cz, ct),
					)
				end
			end

			# Create datasets for GW tendencies.
			for field in (:dudt, :dvdt, :dthetadt)
				if field in output_variables
					create_dataset(
						file,
						string(field),
						datatype(Float32),
						dataspace(
							(x_size, y_size, z_size, 0),
							(x_size, y_size, z_size, -1),
						);
						chunk = (cx, cy, cz, ct),
					)
				end
			end

			if !(typeof(state.namelists.ice.icesetup) <: NoIce) && typeof(state.namelists.ice.cloudcover) <: CloudCoverOn

				# Create datasets for SgsGW variables.
				for field in fieldnames(SgsGW)
					create_dataset(
						file,
						string(field),
						datatype(Float32),
						dataspace(
							(sizex2, sizey2, sizez2, 0),
							(sizex2, sizey2, sizez2, -1),
						);
						chunk = (cx2, cy2, cz2, ct),
					)
				end

				# Create datasets for SgsPredictands
				for field in ("sn", "sq", "sqv")
					create_dataset(
						file,
						field,
						datatype(Float32),
						dataspace(
							(sizex2, sizey2, sizez2, 0),
							(sizex2, sizey2, sizez2, -1),
						);
						chunk = (cx2, cy2, cz2, ct),
					)
				end

				# Create datasets for SgsAuxiliaries
				for field in fieldnames(SgsAuxiliaries)
					create_dataset(
						file,
						string(field),
						datatype(Float32),
						dataspace(
							(sizex2, sizey2, sizez2, 0),
							(sizex2, sizey2, sizez2, -1),
						);
						chunk = (cx2, cy2, cz2, ct),
					)
				end

			end
		end

		# Return.
		return
	end

	# Return.
	return
end
