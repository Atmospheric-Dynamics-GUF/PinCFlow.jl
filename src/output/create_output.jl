function create_output(state::State)

  # Get all necessary fields.
  (; sizex, sizey, sizez) = state.namelists.domain
  (; master) = state.domain

  # Prepare the output on master.
  if master

    # Create the dataset.
    NCDataset("pincflow_output.nc", "c"; share = true) do dataset

      # Define the dimensions.
      defDim(dataset, "x", sizex)
      defDim(dataset, "y", sizey)
      defDim(dataset, "z", sizez)
      defDim(dataset, "t", Inf)

      # Define the dimension variables.
      defVar(dataset, "x", Float32, ("x",))
      defVar(dataset, "y", Float32, ("y",))
      defVar(dataset, "z", Float32, ("x", "y", "z"))
      defVar(dataset, "t", Float32, ("t",))

      # Define the background variables.
      defVar(dataset, "rhobar", Float32, ("x", "y", "z"))
      defVar(dataset, "thetabar", Float32, ("x", "y", "z"))
      defVar(dataset, "n2", Float32, ("x", "y", "z"))
      defVar(dataset, "p", Float32, ("x", "y", "z"))

      # Define the prognostic variables.
      defVar(dataset, "rhop", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "u", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "us", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "v", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "vs", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "w", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "ws", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "wtfc", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "wstfc", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "thetap", Float32, ("x", "y", "z", "t"))
      defVar(dataset, "pip", Float32, ("x", "y", "z", "t"))

      # Return.
      return
    end
  end

  # Return.
  return
end
