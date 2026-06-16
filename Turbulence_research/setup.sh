#!/bin/bash

# Precompile.
julia --project=examples -e 'using Pkg; Pkg.precompile()'
