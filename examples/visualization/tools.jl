using NCDatasets

"""
Composite type for model output.
"""
struct ModelOutput{P <: Dict, A <: Dict, D <: Dict, V <: Dict, G <: Dict}
  parameters::P
  attributes::A
  dimensions::D
  variables::V
  groups::G
end

"""
Import model output data from PinCFlow.
"""
function ModelOutput(path::String)

  # Read namelists.
  lines = readlines(path * "/namelists.txt")

  # Remove white spaces, empty lines, comments and namelists titles.
  lines = [
    strip(split(entry, "!"; limit = 2)[1]) for entry in lines if
    length(strip(entry)) > 1 && !(strip(entry)[1] in ("!", "&"))
  ]

  # Get keys and values.
  keys = [
    lowercase(strip(split(entry, "="; limit = 2)[1])) for
    entry in lines if occursin("=", entry)
  ]
  values = [
    strip(split(entry, "="; limit = 2)[2]) for
    entry in lines if occursin("=", entry)
  ]
  values = [rstrip(entry, ',') for entry in values]

  # Adjust boolen entries.
  values = [replace(entry, r"\bT\b|\.true\." => "true") for entry in values]
  values = [replace(entry, r"\bF\b|\.false\." => "false") for entry in values]

  # Adjust strings.
  values = [replace(entry, "\'" => "\"") for entry in values]

  # Evaluate entries.
  values = [eval(Meta.parse(entry)) for entry in values]

  # Concatenate array elements.
  values = [
    occursin("(", keys[m]) ?
    [
      values[n] for
      n in 1:length(keys) if split(keys[n], "(")[1] == split(keys[m], "(")[1]
    ] : values[m] for m in 1:length(keys)
  ]
  keys = [split(entry, "(")[1] for entry in keys]

  # Save to dictionary.
  parameters = Dict(keys[n] => values[n] for n in 1:length(keys))

  # Import data.
  data_set = NCDataset(path * "/pincflow_data_out.nc")
  attributes =
    Dict(key => data_set.attrib[key] for key in Base.keys(data_set.attrib))
  dimensions = Dict(key => data_set.dim[key] for key in Base.keys(data_set.dim))
  variables = Dict(data_set)
  groups = Dict(
    key =>
      ModelOutput(Dict(), Dict(), Dict(), Dict(data_set.group[key]), Dict())
    for key in Base.keys(data_set.group)
  )

  # Return data.
  return ModelOutput(parameters, attributes, dimensions, variables, groups)
end
