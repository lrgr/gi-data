# Genetic Interaction datasets

Genetic Interaction datasets for various species. Datasets are processed taking into account various idiosyncracies of protocols and studies (essential vs. non-essential genes, P-values, TS-alleles etc.).

## 1. Usage

We recommend users download data from static links provided in the provided YAML in `/example`. (Datasets can be downloaded as in `/example`. Unless specified otherwise, data is serialized and stored with cloudpickle (See example for usage).

For maintainers and advanced user, see 3. to build data from source.

## 2. General information
Processed GI datasets:
- Collins et al. Nature 2007
- Roguev et al. Science 2008
- Costanzo et al. Science 2010
- Costanzo et al. Science 2016
- Ryan et al. Molecular Cell 2012

## Build data from source (i.e. maintainers)
1. First build each dataset in respective directories
2. Run `snakemake all` to collect all datasets for release/upload
3. Run `snakemake verify` to verify built files against committed md5 checksums
4. (For maintainers only: run `sh umdobj_upload.sh <ver>` to upload to UMIACS ObjectStore)

**Note:** a local python package is installed (see environment.yml) for convenience and standardization of dataset serialization.
