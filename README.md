# Summary

Plots and analysis for [Zhou _et al._, _Noncoding mutations target cis-regulatory elements of the FOXA1 plexus in prostate cancer_, 2019]().

Data used for these analyses can be found in `data/`.
Analyses performed can be found in `results/`.

## Setup

For reproducibility, the core packages (and their versions) used for these analyses can be found in `environment.yaml`.
To reproduce this environment, use [Conda](https://conda.io/docs/) with the following command:

```shell
conda env create --name Zhou2019 --file environment.yaml
```

This repository can be downloaded using `git clone`.

## Analysis

To re-run the analysis in a specific folder (i.e. `results/<analysis>/`), simply run `snakemake` from within this folder, or follow the `README.md` files.
