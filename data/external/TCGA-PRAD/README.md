# Summary

This folder contains RNA-seq FPKM values from TCGA for the prostate cancer project.
I've downloaded these from the GDC Data Portal with the following filters:

| Parameter     | Filter                         |
| ------------- | ------------------------------ |
| Project Id    | TCGA-PRAD                      |
| Data Category | Transcriptome Profiling        |
| Data Type     | Gene Expression Quantification |
| Workflow Type | HTSeq - FPKM                   |

The manifest for this dataset is `gdc_manifest_20180717_162017.txt`.

I've downloaded the data using the `gdc-client` tool, and the manifest within this folder.
This was done via the rule defined in `Snakefile` and the `snakemake` command.

## Notes

All of the FPKM text files have the same 2-column format: `ENSEMBL Gene ID` and `FPKM`