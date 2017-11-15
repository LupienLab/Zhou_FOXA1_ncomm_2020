# Summary

Running SMURF on prostate cancer mutations.
Mutations come from 200 intermediate-risk primary prostate cancer samples.

Regions of interest are DHS regions from LNCaP cells, intersected with catalogue of H3K27ac peaks from 20 patients (from the above cohort).

## Steps

1. Intersect LNCaP DHS regions with H3K27ac catalogue
    ```shell
    bedtools intersect -a ../../data/processed/LNCaP_DHS.bed -b ../../data/processed/20PCa_H3K27acpeaks_extended500bp_merged.bed -sorted > regions_LNCaP-DHS_H3K27ac-catalogue.bed
    ```
1. Run SMURF
    ```shell
    sh run-smurf.sh
    ```