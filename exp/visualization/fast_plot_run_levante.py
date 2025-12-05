#!/usr/bin/env python3
#python3 test/fast_plot_run_levante.py
"""
Python rewrite of test/fast_plot_run_levante.jl
Reads `mountain_wave.h5`, plots a 2D pcolormesh and a PDF of the selected field.
"""
import os
from pathlib import Path
import h5py
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi

# Base font-size (points). Change this to scale all plot text.
FS = 12

INPUT = Path.home() /"spr/ice_mountain_wave.h5"
OUT_DIR = Path.home() / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

with h5py.File(INPUT, "r") as data:
    # Read grid and fields (units converted to km)
    x_vals = np.array(data["x"]) * 0.001
    y_vals = np.array(data["y"]) * 0.001
    z = np.array(data["z"]) * 0.001  # expected shape (nx, ny, nz)

    n_y = data["w"].shape[2]
    iy = 0 #n_y // 2
    iys = 0
    print("iy =", iy)
    print(z[0,0,:])

    tidx = 20
    fld = np.array(data["u"][tidx, ...])

    print("size fld ", fld.shape)
    print("fld ", np.nanmax(fld), np.nanmin(fld))

    # --- Plot 2D field (x vs z) at slice iy ---
    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(2, 1, 1)

    # Slice arrays to shapes (nx, nz) for plotting
    F = fld[:, iy, :]

    #pcm = ax1.pcolormesh(X, Z, Z, cmap="Blues")
    #pcm = ax1.pcolormesh(x_vals, z[:,iy,:], F, cmap="Blues")
    pcm = ax1.contourf(np.broadcast_to(x_vals[None,:], F.shape), z[:,iy,:], F, levels=10)
    cbar = fig.colorbar(pcm, ax=ax1)
    #cbar.ax.tick_params(labelsize=FS)
    ax1.set_title("data: Res", fontsize=1.5 * FS)
    ax1.set_xlabel("x (km)", fontsize=1.2 * FS)
    ax1.set_ylabel("z (km)", fontsize=1.2 * FS)
    ax1.tick_params(labelsize=FS)

    # The second subplot for the secondary dataset was disabled in the Julia file.
    # We'll leave the second subplot blank to match the structure.
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.axis("off")

    plt.tight_layout()
    fig.savefig(OUT_DIR / "mountain_wave.png")
    plt.close(fig)
