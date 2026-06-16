#!/bin/bash

set -e

julia --project=Turbulence_research -e 'include("Turbulence_research/Turbulence_research.jl"); using .Turbulence_research; wkb_wp_1d(; x_size=1, y_size=1, z_size=200, tmax = 2.0, output_interval = 1.0)' &> wkb_wp_1d.log