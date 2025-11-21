dirname=$1
numberprc=$2

mpiexec=$(julia --project=. -e 'using MPICH_jll; println(MPICH_jll.mpiexec_path)')
${mpiexec} -n ${numberprc} julia --project test/ice_${dirname}.jl &&
mkdir -p ../runs/${dirname} &&
mv test/pincflow_output.h5 ../runs/${dirname} && cp -pr src test/ice_${dirname}.jl ../runs/${dirname}/.