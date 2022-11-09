[![DOI](https://zenodo.org/badge/409672778.svg)](https://zenodo.org/badge/latestdoi/409672778)

# Quality control and removal of technical variation of NMR metabolic biomarker data in ~120,000 UK Biobank participants

This repository houses code used to identify and remove technical variation form the NMR metabolic biomarker data currently (2021) available in ~120,000 UK Biobank participants. For details, please see our preprint, Ritchie S. C. *et al.*, Quality control and removal of technical variation of NMR metabolic biomarker data in ~120,000 UK Biobank participants, **medRxiv** (2021). doi: [10.1101/2021.09.24.21264079](https://www.medrxiv.org/content/10.1101/2021.09.24.21264079v1).

## Source code organisation

Source code is located in the `src/` folder and organised into sequential steps and utility functions:

### 1. Sample quality control

Scripts located under `src/01_sample_harmonization/` perform sample quality control of the pre-release raw data made available to early access analysts under UK Biobank project 30418 and subsequent harmonization of sample QC with the samples made available by UK Biobank in their public release.

- `01_raw_nightingale_sample_qc.R`: performs sample quality control of the pre-release raw data.
- `02_extract_UKB_nmr.sh`: extracts the raw data from decoded UK Biobank download of the publicly released raw data.
- `03_harmonize_sample_QC.R`: identifies the samples passing QC in the pre-release raw data that were also made available by UK Biobank in the public download.

### 2. Biomarker quality control and diagnostic plots

Scripts located under `src/02_biomarker_qc/` do the following:

- `01_technical_qc.R`: Removes the effects of technical variation from biomarkers and generates step-wise intermediate data at each step of the process.
- `02_diagnostic_plots.R`: Generates diagnostics plots for all biomarkers showing their relationship with technical covariates at each step in the removal of technical variation and how concentrations change at each step.
- `03_qc_no_rederivation.R`: Generates data for comparing the main approach to one in which composite biomarkers and ratios are directly adjusted instead of recomputed after adjustment.
- `04_age_sex_bmi.R`: Taking the data with technical variation removed (output by `01_technical_qc.R`) creates and additional dataset adjusted for age, sex, and BMI, both with and without recomputing composite biomarkers and ratios, so that the effects of direct adjustment vs. recomputation of composite biomarkers and ratios after adjustment can be more effectively compared.

### 3. Scripts to generate paper figures and tables

Scripts located under `src/03_paper/` generate figures and tables for the paper, including downstream analyses (e.g. example GWAS). Scripts are numbered sequentially, and generally correspond to the order of figures and tables (main and supp) as presented in the paper. Some figures may have components from multiple scripts.

### 4. Scripts to check ukbnmr R package and new UK Biobank releases

The script `src/04_pkg_sanity_check/01_pkg_sanity_check.R` sanity checks the results of using the [ukbnmr](https://github.com/sritchie73/ukbnmr) R package to adjust raw data released from UK Biobank for technical variation, particularly as the data released by UK Biobank differs in sample content from the pre-release data: blind duplicate samples are not included, and sample inclusion may change over time due to participant withdrawals.

#### Utility functions

Scripts stored in `src/utilities/` contain general purpose functions which may be called from other scripts above. Currently the only script is `logit.R` which contains functions for performing logit transformation (the analog of log transformation for percentages) and for inversing the function (so that logit distributions may be converted back to % units).

## Dependencies

### Software dependencies

Here we list the software and versions used throughout the pipeline. You will need the same software and packages to use these scripts. We recommend also that you use at least the same version numbers or higher, as code has not been tested with earlier versions.

#### 1. Sample quality control

- R version 4.0.3 (2020-10-10)
- R packages:
  - data.table version 1.13.2
  - ukbnmr version 0.3.0 # Currently installable via `remotes::install_github("sritchie73/ukbnmr", ref="development")`
- UK Biobank's ukbconv tool (download from http://biobank.ndph.ox.ac.uk/ukb/download.cgi) stored in `src/ukbtools/`.

#### 2. Biomarker quality control and diagnostic plots

- R version 4.0.3 (2020-10-10)
- R packages:
  - data.table version 1.13.2
  - ukbnmr version 0.3.0 # Currently installable via `remotes::install_github("sritchie73/ukbnmr", ref="development")`
  - MASS version 7.3-53
  - ggplot2 version 3.3.2
  - ggthemes version 4.2.4
  - palettetown version 0.1.1
  - readstata13 version 0.9.2

The `readstata13` package is optional. It is used in `04_age_sex_bmi.R` to load in the UK Biobank phenotype dataset curated by the Cardiovascular Epidemiology Unit at the University of Cambridge, which is stored in STATA 13 format and contains the age, sex, and BMI information for participants in project 7439. You will need to replace this with your own extracted dataset and readstata13 can be omitted unless you are also working with a STATA 13 file.

#### 3. Scripts to generate paper figures and tables

- R version 4.0.3 (2020-10-10)
- R packages:
  - data.table version 1.13.2
  - ukbnmr version 0.3.0 # Currently installable via `remotes::install_github("sritchie73/ukbnmr", ref="development")`
  - MASS version 7.3-53
  - ggplot2 version 3.3.2
  - ggrastr version 0.2.3
  - ggnewscale version 0.4.3
  - cowplot version 1.1.0
  - ggthemes version 4.2.4
  - palettetown version 0.1.1
  - RColorBrewer version 1.1-2
  - rcartocolor version 2.0.0
  - palr version 0.2.0
  - ochRe version 1.0.0 # GitHub package only, install via `remotes::install_github("ropenscilabs/ochRe")`
  - WGCNA version 1.69 # Install from BioConductor, `BiocManager::install("WGCNA")`
  - foreach version 1.5.1
  - doMC version 1.3.7
  - httr version 1.4.2
  - jsonlite version 1.7.1
  - xml2 version 1.3.2  
  - annotables version 0.1.91 # Install from BioConductor, `BiocManager::install("annotables")`
  - readstata13 version 0.9.2
  - survival version 3.2-7
- Plink 2 version v2.00a3LM AVX2 Intel (2 Mar 2021) # downloadable from https://www.cog-genomics.org/plink/2.0/

The `readstata13` package is optional. It is used in `14_gwas.R`, `15_cad_stroke_episcan.R`, and `17_sample_pca.R` to load in the UK Biobank phenotype dataset curated by the Cardiovascular Epidemiology Unit at the University of Cambridge, which is stored in STATA 13 format. This dataset is loaded to obtain participant age, sex, BMI, ethnicity, fasting time, and lipid lowering medication usage.You will need to replace this with your own extracted dataset and readstata13 can be omitted unless you are also working with a STATA 13 file.

The scripts `01_tech_covariate_r2.R`, `14_gwas.R`, and `15_cad_stroke_episcan.R` are parallized. The parallel backend is registered using the `doMC` package. If you are using a windows system, this package will not work for you, and you will need instead to modify the code to register a parallel backend using the `doParallel` package instead. You will also need to modify the number of cores to be compatible with your system. `14_gwas.R` also detects the number of cores automatically assuming the script was submitted as a job on a SLURM job submissions system. You may need to instead hard code the number of desired cores.

The script `14_gwas.R` also extensively makes use of the `system()` function to run commands in the BASH shell environment. These will not work on a windows system, but should be portable across linux distros and OSX. Those using Windows 10 should be able to execute this in the Windows Subsytem for Linux environment. The BASH version we used was GNU bash, version 4.2.46(2)-release (x86\_64-redhat-linux-gnu) on RedHat Scientific Linux release 7.9 (Nitrogen).

### Paper drafting and figure annotation

Microsoft Office Professional Plus 2016 was used to draft the manuscript (Microsoft Word) and curate supplemental tables (Microsoft Excel) on Windows 10 Enterprise edition. Inkscape version 0.92.3 was used to layout and annotate figures from the figure components generated within the R scripts. 

### Data folders

The scripts contain a variety of hard coded links to files stored under the `data/` directory. This folder is not provided in this GitHub repository as we cannot make participant level data publicly available. You will need to create these yourself with data obtained from UK Biobank directly. See https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access for further details.

- `data/raw/nightingale/`: folder containing pre-release raw data made available to early access researchers under project 30418. These are not available from the UK Biobank showcase, and will need to be specifically requested in addition to the raw NMR metabolic data provided through UK Biobank. Note the file `Cambridge_Nightingale_bridge.tsv` contains a mapping between internal sample IDs (not project specific) and project 7439 participant IDs. You will also need to obtain a mapping file from UK Biobank to your relevant project here.
- `data/raw/ukbiobank/decoded/`: folder containing raw UK Biobank data that has been decoded using your project's encryption key. Files in this folder are given project specific alphanumeric codes, which you will need to change in the relevant scripts.

For downstream analyses (GWAS and association scans for incident stroke and coronary artery disease) and comparison to clinical biochemistry a variety of predefined files from other projects have been used:

- `data/imputed_pgen/`: Imputed genotype data from UK Biobank converted to plink2 binary format.
- `data/genetic_reference/`: reference files for UK Biobank genetic data, including `ukb_sqc_v2.txt` and files output by `qctool2` summarising variant information. 
- `data/Berisa_2016_LD_blocks/EUR_1000G_ind_ld_blocks.bed`: Pre-defined LD blocks for european ancestry from Berisa *et al.* 2016 obtained from https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier\_ls-all.bed
- `data/ceu_curated_phenotypes/`: Department curated UK Biobank phenotype data in STATA13 format (see above).
- `data/my_curated_phenotypes/biomarkers/output/biomarkers_fixed_limits.txt`: File containing clinical biochemistry data from UK Biobank with values above and below detection limit extracted from the respective fields and replace with minimum and maximum observed values in the data.
- `data/my_curated_phenotypes/endpoints/`: Folder containing time-to-first-event data for endpoints defined from ICD10 codes and prevalent disease defined from ICD10, ICD9, and self-reported history.


