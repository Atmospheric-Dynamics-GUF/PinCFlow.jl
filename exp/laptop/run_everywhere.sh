#!/usr/bin/env bash
set -euo pipefail

# run_everywhere.sh - determine machine locations and prepare run directories
# Usage: run_everywhere.sh [runName]

runName=${1:-default_run}

# Determine host and user
host_name=$(hostname)
user_name=$(whoami)

if [[ "$host_name" == *levante* ]]; then
  # Levante cluster
  data_path="/scratch/b/${user_name}/PF/runs"
  dirHome="/home/${user_name}/PF/pinc"
  dirScratch="/scratch/b/${user_name}/PF/runs/${runName}"
  dirSaveCode="/home/${user_name}/PF/code_runs/${runName}"
  #reference_path="/scratch/b/${user_name}/PF/pinc/reference"
elif [[ "$host_name" == *login* ]]; then
  # Goethe cluster (login node)
  data_path="/scratch/atmodynamics/${user_name}/PF/runs"
  # Define commonly used directories (adjust user/home base if needed)
  dirHome="/home/atmodynamics/${user_name}/PF/pinc"
  dirScratch="/scratch/atmodynamics/${user_name}/PF/runs/${runName}"
  dirSaveCode="/home/atmodynamics/${user_name}/PF/code_runs/${runName}"
  #reference_path="/scratch/atmodynamics/${user_name}/PF/pinc/reference"
else
  # Local machine or unknown
  dirHome="/home/dolaptch/PF/pinc"
  dirScratch="/home/dolaptch/PF/runs/${runName}"
  dirSaveCode="/home/dolaptch/PF/code_runs/${runName}"
  #reference_path="../"
fi

mkdir -p "${dirScratch}"
mkdir -p "${dirSaveCode}"

namelist=${dirHome}/exp/namelists_run.jl

# Copy source tree and this script to the save directory
if [[ -d "${dirHome}" ]]; then
  cp -rp "${dirHome}/src" "${dirHome}/Project.toml" "${dirHome}/Manifest.toml" "${namelist}" "${dirSaveCode}/."
fi
cp -p "${BASH_SOURCE[0]}" "${dirSaveCode}/." || true

# Clean the scratch working directory (CAREFUL: this removes files)
if [[ -d "${dirScratch}" ]]; then
  # remove contents only
  rm -rf "${dirScratch:?}"/* || true
fi

#cd "${dirHome}"
cd "${dirScratch}"

# Placeholder: execute the Julia script or other commands here
# Example: julia --project=${dirHome} exp/fast_plot_run.jl
julia --project="${dirHome}" "${namelist}"

exit 0