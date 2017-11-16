# Summary

Running SMURF on prostate cancer mutations.
Mutations come from 200 intermediate-risk primary prostate cancer samples.

Regions of interest are DHS regions from LNCaP cells, intersected with catalogue of H3K27ac peaks from 20 patients (from the above cohort).

This folder uses data that was used to generate the results in `/results/1A_smurf/`, but there were some questions about how that data was generated.
So this folder is to explicitly regenerate the data and those figures.
These results should be used instead of the figures from `/results/1A_smurf/`.

## Steps

1. Intersect LNCaP DHS regions with H3K27ac catalogue
    ```shell
    # intersect BED files
    bedtools intersect -a ../../data/processed/LNCaP_DHS.bed -b ../../data/processed/20PCa_H3K27acpeaks_extended500bp_merged.bed -sorted > regions_LNCaP-DHS_H3K27ac-catalogue.temp.bed
    # add 4th column required for SMURF
    awk -v FS="\t" -v OFS="\t" '{print $1, $2, $3, "."}' regions_LNCaP-DHS_H3K27ac-catalogue.temp.bed > regions_LNCaP-DHS_H3K27ac-catalogue.bed
    rm regions_LNCaP-DHS_H3K27ac-catalogue.temp.bed
    ```
1. Run SMURF
    ```shell
    sh run-smurf.sh
    ```