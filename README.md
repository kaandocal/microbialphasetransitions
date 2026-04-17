This repository contains Julia code for the paper:

K. Öcal, S. Ghrabli, M. P. H. Stumpf, "Phase transitions in microbial lineage trees" (2025)

The code is based on [this repository](https://github.com/kaandocal/twoclock). It consists of `sim.jl`, which is used to simulate populations with a fixed carrying capacity, and a model `plasmids.jl`.

To generate data, run the file `run.jl` as follows:
```
julia run.jl --tmax [TMAX] --N [N] [α] [β] 
```

This will generate output in the form of a JLD2 file in the subfolder "data". The files `plot_plasmid.jl` and `plot_heatmap.jl` can be used to plot the growth rates and selection coefficients as in the paper -- this requires the data to be simulated first.

If you have any questions or comments, please feel free to contact the authors or open a pull request on GitHub.
