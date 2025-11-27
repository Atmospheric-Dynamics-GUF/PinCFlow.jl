dirname=$1
option=$2
#julia ${option} --project exp/ice_${dirname}.jl &&
mkdir -p ../runs/${dirname} &&
mv exp/pincflow_output.h5 ../runs/${dirname} && cp -pr src exp/ice_${dirname}.jl ../runs/${dirname}/.