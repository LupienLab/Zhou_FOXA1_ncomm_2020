# Summary

This folder contains data for a similar test of essentiality of FOXA1, but in luminal breast cancer as opposed to prostate cancer.

The permutation tests are run in the same manner as those in `1C_gene-essentiality/`, just using breast cancer cell lines based on their subtype classification according to work from Wail.

## Results

### "Very Good" classification

These are cell lines that are classified consistently across 3 breast cancer subtyping algorithms: PAM50, SCMOD2, and SCMGENE.

There are 3 cell lines that meet this criteria for the luminal subtype:

* T-47D
* MCF-7
* EFM-19

T-47D has the lowest FOXA1 essentiality score (i.e. most essential) across all cell lines, meaning the results from the CRISPRi data will give a p-value of 0.

For the RNAi data, selecting a single cell line at random (Test 1) less than EFM-19 or MCF-7 is 2.78%, and 2 cell lines at random is 13.65% (Test 2).
These are similar numbers to those found for FOXA1 in prostate cancer cell lines.

### "Good" classification

These are cell lines that are classified consistently between PAM50, SCMOD2 (this includes the "very good" cell lines).
Unfortunately, there are no more luminal cell lines with "good" classification, so the result will be the same as above.

### All breast cancer cell lines

There are 2 breast cancer lines in the CRISPRi data:

* CAL-120
* T-47D

There are also 13 breast cancer lines in the RNAi data:

* BT-20
* BT-474
* CAL-120
* CAL-51
* EFM-19
* HCC-1187
* HCC-1395
* HCC-1954
* HCC-2218
* HCC-70
* MCF-7
* MDA-MB-453
* ZR-75-30

Test 1 is possible for both datasets, but Test 2 is intractable for the RNAi data (> 2.4e20 combinations).

For the CRISPR data, the probability of selecting a single cell line with FOXA1 essentiality lower than one of the breast cancer lines is 43.9%, and the probability of selecting 2 with a lower median essentiality than that of the breast cancer lines is 7.58%.

For the RNAi data, the probability of selecting a single cell line with FOXA1 essentiality lower than one of the breast cancer lines is 26.8%.

## Conclusions

FOXA1 is particularly essential in luminal breast cancer, much like prostate cancer, but this is not the case across all breast cancers, more generally.