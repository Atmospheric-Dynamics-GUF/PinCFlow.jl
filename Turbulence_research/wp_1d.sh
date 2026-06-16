#!/bin/bash

set -e

julia --project=Turbulence_research -e 'include("Turbulence_research/Turbulence_research.jl"); using .Turbulence_research; wp_1d(; x_size=16, y_size=16, z_size=400, tmax = 2.0, output_interval = 1.0)' &> wp_1d.log