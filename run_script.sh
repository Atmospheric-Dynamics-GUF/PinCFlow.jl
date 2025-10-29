dirname=$1
julia --project test/ice_$dirname.jl &&
mkdir -p ../runs/$dirname &&
mv test/pincflow_output.h5 ../runs/$dirname && cp -pr src ../runs/$dirname/.
    