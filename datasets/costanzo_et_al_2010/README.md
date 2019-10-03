
# Costanzo et al. 2010 data

Website for supplement: http://boonelab.ccbr.utoronto.ca/supplement/costanzo2009/

## Usage:
Run `snakemake all` to produce GI files for non-essential genes. Different quantitative measurements are saved to `processed` with sensible names.

- `costanzo2010-GI_scores.pkl` contains GI scores
- `costanzo2010-pval.pkl` contains P-values
- Other files contain single mutant, and wild-type fitnesses

## Details:
- For Non-essential genes/deletion mutants only:
    - Sparsity: 0.09863028294400844
    - Shape: (1377, 3885)
- From the logs...
    - Raw GIs shape (including DAMP, TS alleles): (1711, 3885)
        - Removing DAMP genes
                - 0 DAMP genes on cols
                - 120 DAMP genes on rows
        - GIs shape: (1591, 3885)
        - Removing TS genes
                - 0 TS genes on cols
                - 214 TS genes on rows
        - GIs shape: (1377, 3885)
        - GIs shape: (1377, 3885)
        - Sparsity: 0.09863028294400844
        - Interactions measured : 4822008
        - Interactions missing: 527637

## Kinds of released GI files:

Different 'cutoffs' for epsilon scores and p-values were selected to count the number of 'actual' GIs.

### Raw data:
All GI scores

### Lenient Cuttoff:
The file contains the SGA genetic interaction dataset with a lenient cutoff applied (p-value < 0.05). Reciprocal interactions (AB vs BA) were processed as follows: if AB and BA show opposite interaction signs (AB is positive and BA is negative, or viceversa), both pairs were removed; if AB and BA show the same interaction sign (both positive or both negative), the interaction with the lowest p-value was retained and both pairs are reported with that interaction. The file is provided in a tab-delimited format with 7 columns:

- ~700K pairs

### Intermediate Cuttoff:
The file contains the SGA genetic interaction dataset with an intermediate cutoff applied (|ε| > 0.08, p-value < 0.05). Reciprocal interactions (AB vs BA) were processed as follows: if AB and BA show opposite interaction signs (AB is positive and BA is negative, or viceversa), both pairs were removed; if AB and BA show the same interaction sign (both positive or both negative), the interaction with the lowest p-value was retained and both pairs are reported with that interaction. The file is provided in a tab-delimited format with 7 columns:

### Stringent Cuttoff:
The file contains the SGA genetic interaction dataset with a stringent cutoff applied (ε < -0.12, p-value < 0.05 or ε > 0.16, p-value < 0.05). Reciprocal interactions (AB vs BA) were processed as follows: if AB and BA show opposite interaction signs (AB is positive and BA is negative, or viceversa), both pairs were removed; if AB and BA show the same interaction sign (both positive or both negative), the interaction with the lowest p-value was retained and both pairs are reported with that interaction. The file is provided in a tab-delimited format with 7 columns:

## Other details:
- Biological process annotations are [here](http://boonelab.ccbr.utoronto.ca/supplement/costanzo2009/bioprocess_annotations_costanzo2009.xls)
- Temperature sensitive genes and DAMP-ed genes are marked with `_DAMP` and `_tsq{n}` suffixes.