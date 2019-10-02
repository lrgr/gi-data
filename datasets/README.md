# Procesesd GI datasets

## Datasets:
- Collins et al. Nature 2007
- Roguev et al. Science 2008

## Usage:
1. First build each dataset in respective directories
2. Run `snakemake all` to collect all datasets for release/upload
3. Run `snakemake verify` to verify built files against committed md5 checksums
4. (For maintainers only: run `sh umdobj_upload.sh <ver>` to upload to UMIACS ObjectStore)