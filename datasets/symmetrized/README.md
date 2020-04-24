# Associating unique gene pairs to unique GIs.

Gene pairs in (Array, Query) position, (A,B) or (B,A) can have different GI scores.

To assign unique gene pairs to unique GI scores we process GI scores as follows.
- If (A,B) and (B,A) have different signs, discard GI scores
- Else, if scores have P-values, retain score with lower P-value
- Else, retain the mean score.

## Usage:

To process Costanzo et al., (Science 2010), and Ryan et al. (Molecular Cell, 2012), run:

    snakemake all

