# Data from Ryan CJ et al. Molecular Cell, 2012


## Details
From the paper:
- Screened 953 alleles (Table S1 available online) of 876 genes against a fission yeast mutant library containing more than 2,000 deletions (Table S1)
- resulting in an E-MAP containing ∼1.6 million pairwise measurements (Data Sets S1 and S2). 
- The majority of the genes screened are broadly conserved across eukaryotes, with subsets that are fungal- and fission yeast-specific (Figure 1A and Table S1).

## Parsed data
Parsed data in `output/processed`:
- `ryan_2012_sc_merged_gis.cpkl`
    - DAMP and deletion alleles for Sc (I think...)
    - 52% sparse
    - 4439 x 4439 unique genes. The dataset is square and symmetric
    - For Sc data, there is one error. The entry `YLR176`, is misnamed, and is fixed to be `YLR176C`

- `ryan_2012_sp_gis.cpkl`
    - Sp EMAP (may have damped alleles...) 
    - 16% sparse
    - 862 x 1955 unique genes. (407 genes in both query and array)

## Issues:
- In S.c data, aliased names of the same gene appears in different rows/cols in the GI matrix. Not sure why this is. These aliases are found heuristically by checking if genes in rows are prefixed, by name, by genes in the columns, and vice versa. (These aliased gene pairs are printed to the log...)