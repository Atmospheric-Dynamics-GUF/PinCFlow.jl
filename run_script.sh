dirname=$1
option=$2
julia ${option} --project test/ice_${dirname}.jl &&
mkdir -p ../runs/${dirname} &&
mv test/pincflow_output.h5 ../runs/${dirname} && cp -pr src ../runs/${dirname}/.