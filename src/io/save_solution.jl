function save_solution_file_2d(grid_array,
                               vars,
                               varnames;
                               dt = 0.0,
                               t = 0.0,
                               equation_name = "Unspecified Equations",)
    nx, nz = size(grid_array)

    h5open(filename, "w") do file
        # Add context information as attributes
        attributes(file)["ndims"] = 2
        attributes(file)["equations"] = "Unspecified equations"
        attributes(file)["polydeg"] = 0
        attributes(file)["n_vars"] = n_vars
        attributes(file)["n_elements"] = nx * ny
        attributes(file)["mesh_type"] = "StructuredMesh" # For Trixi2Vtk
        attributes(file)["mesh_file"] = "mesh.h5"
        attributes(file)["time"] = convert(Float64, time) # Ensure that `time` is written as a double precision scalar
        attributes(file)["dt"] = convert(Float64, dt) # Ensure that `dt` is written as a double precision scalar
        attributes(file)["timestep"] = iter

        # Store each variable of the solution data
        var_names = ("Density", "Velocity x", "Velocity y", "Pressure")
        for v in 1:n_vars
            # Convert to 1D array
            file["variables_$v"] = vec(data[v, .., :])

            # Add variable name as attribute
            var = file["variables_$v"]
            attributes(var)["name"] = var_names[v]
        end

        # Store element variables
        for (v, (key, element_variable)) in enumerate(element_variables)
            # Add to file
            file["element_variables_$v"] = element_variable

            # Add variable name as attribute
            var = file["element_variables_$v"]
            attributes(var)["name"] = string(key)
        end
    end
end
