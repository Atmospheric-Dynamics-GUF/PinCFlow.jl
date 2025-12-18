#!/usr/bin/env bash
set -euo pipefail

# run_everywhere.sh - determine machine locations and prepare run directories
# Usage: run_everywhere.sh [runName]

#parallel="yes"
numberprc=2
optimization="-O3"

parallel="${3:-no}"
runName=${1:-default_run}
namelist_run=${2:-namelists_run};

# Determine host and user
host_name=$(hostname)
user_name=$(whoami)

if [[ "$host_name" == *levante* ]]; then
  # Levante cluster
  dirHome="/home/b/${user_name}/PF/pinc"
  dirScratch="/scratch/b/${user_name}/PF/runs/${runName}"
  dirSaveCode="/home/b/${user_name}/PF/code_runs/${runName}"
elif [[ "$host_name" == *login* ]]; then
  # Goethe cluster (login node)
  dirHome="/home/atmodynamics/${user_name}/PF/pinc"
  dirScratch="/scratch/atmodynamics/${user_name}/PF/runs/${runName}"
  dirSaveCode="/home/atmodynamics/${user_name}/PF/code_runs/${runName}"
else
  # Local machine or unknown
  dirHome="/home/dolaptch/PF/pinc"
  dirScratch="/home/dolaptch/PF/runs/${runName}"
  dirSaveCode="/home/dolaptch/PF/code_runs/${runName}"
fi

# Clean the scratch working directory 
if [[ -d "${dirScratch}" ]]; then
  # remove contents only
  rm -rf "${dirScratch:?}"/* || true
fi
if [[ -d "${dirSaveCode}" ]]; then
  # remove contents only
  rm -rf "${dirSaveCode:?}"/* || true
fi

mkdir -p "${dirScratch}"
mkdir -p "${dirSaveCode}"

namelist=${dirHome}/exp/${namelist_run}.jl
    
# Copy source and this script to the save directory
if [[ -d "${dirHome}" ]]; then  
  cp -rp "${dirHome}/src" "${dirHome}/Project.toml" "${dirHome}/Manifest.toml" "${namelist}" "${dirSaveCode}/."
fi
cp -p "${BASH_SOURCE[0]}" "${dirSaveCode}/." || true

#cd "${dirHome}"
cd "${dirScratch}"

# Placeholder: execute the Julia script or other commands here
# Example: julia --project=${dirHome} exp/fast_plot_run.jl
if [[ $user_name == "dolaptch" ]]; then
  if [[ ${parallel} == "yes" ]]; then
    mpiexec=$(julia --project="${dirHome}" -e 'using MPICH_jll; println(MPICH_jll.mpiexec_path)')
    ${mpiexec} -n $numberprc julia "${optimization}" --project="${dirHome}" "${namelist}" $numberprc
    exit 0 
  else
    julia "${optimization}" --project="${dirHome}" "${namelist}"
  fi
fi

if [[ $user_name == "b381734" ]]; then
   echo $dirSaveCode
   echo $dirScratch
   cp -p "${dirHome}/exp/levante/ice_dump.sh" "${dirSaveCode}/."
   sbatch ${dirHome}/exp/levante/ice_dump.sh ${namelist}
fi 

exit 0