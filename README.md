Developing serum proteomics based prediction models of disease progression in ADPKD
================
2025-06-05

- [Scripts](#scripts)
- [First Run](#first-run)
- [Building Predictive Models](#building-predictive-models)
- [Validation](#validation)
- [IgAN cohort](#igan-cohort)
- [Enrichment Analysis on Whole Proteome of
  SC](#enrichment-analysis-on-whole-proteome-of-sc)

## Scripts

| Scripts                             | Details                                                                         |
|-------------------------------------|---------------------------------------------------------------------------------|
| 00functions.R                       | all functions                                                                   |
| SC/SC_01_annot.half.R               | clinical annotation                                                             |
| SC/SC_02_interventions.R            | intervention information                                                        |
| SC/SC_03_mergedCreasUpd.R           | merging clinical annotation with intervention information and calculating slope |
| SC/SC_04_annots_merging.R           | merging proteome and clinical annotation                                        |
| SC/SC_05_prots.R                    | proteome data                                                                   |
| SC_10_mainanalysis.R                | main analysis                                                                   |
| SC_20_enrichmentGraph.R             | enrichment analysis for the whole proteome                                      |
| VC/VC_01_ECclinical.R               | clinical annotation for external cohort (EC)                                    |
| VC/VC_02_mergedClinical.Rannotation | merging proteome and clinical                                                   |
| VC/VC_03_proteome_load.R            | loading proteome data                                                           |
| VC/VC_03_proteome_open.R            | proteome data                                                                   |
| VC_10_mainanalysis.R                | validation                                                                      |
| VC_11_furtheranalysis.R             | further analysis on validation                                                  |
| VC_12_IgAN.R                        | IgAN cohort                                                                     |

## First Run

First, start with
“~/Data/slope-models-adpkd/code/VC/VC_02_mergedClinical.R”. This script
will generate a clinical annotation file.

``` r
source("~/Data/slope-models-adpkd/code/VC/VC_02_mergedClinical.R", echo = FALSE)
```

Second, open this script
(“slope-models-adpkd/code/VC/VC_03_proteome_load.R”) and comment out
line 34. Then source the same script. By doing so, an object named
“data_s_gronDat” will be saved into /data folder with the same name as
the object.

``` r
source("~/Data/slope-models-adpkd/code/VC/VC_03_proteome_load.R", echo = FALSE)
```

After saving the object, quit

``` r
quit()
```

## Building Predictive Models

First start with “~/Data/slope-models-adpkd/code/SC_10_mainanalysis.R”.
This script will save the environment and generate figures and tables.

``` r
source("~/Data/slope-models-adpkd/code/SC_10_mainanalysis.R", echo = FALSE)
```

After sourcing, following

- tables are not saved, but can be found in the viewer

  - Table S6, S7, S9, S11, and

- figures and tables are saved in /data/results folder:

  - Table 2, S3, S4, S10

  - Figure 1, 2, S5

  - Figure S10 (from ~/Data/slope-models-adpkd/code/SC/SC_05_prots.R)

  - Data S3, S4, S5



Then, quit

``` r
quit()
```

## Validation

### Validating the Models

First start with “~/Data/slope-models-adpkd/code/VC_10_mainanalysis.R”.
This script will save the environment and generate figures and tables.

``` r
source("~/Data/slope-models-adpkd/code/VC_10_mainanalysis.R", echo = FALSE)
```

After sourcing, following

- figures and tables are saved in /data/results folder

  - Table 1, Table S2 (as Table1.xlsx)

  - Table S8 (as “patientsizeValidation” and “sampsizeValidation” txt
    files)

  - Figure 3, S6, S7, S8

  - Figure S11 (from
    ~/Data/slope-models-adpkd/code/VC/VC_03_proteome_load.R)

Then, quit

``` r
quit()
```

### Further Analyses on Validation

First start with
“~/Data/slope-models-adpkd/code/VC_11_furtheranalysis.R”. This script
will generate a figure and a table.

``` r
source("~/Data/slope-models-adpkd/code/VC_11_furtheranalysis.R", echo = FALSE)
```

After sourcing, following figure and table are saved in /data/results
folder

- Table S1

- Figure S9

Then, quit

``` r
quit()
```

## IgAN cohort

First start with “~/Data/slope-models-adpkd/code/VC_12_IgAN.R”. This
script will generate a figure and a table.

``` r
source("~/Data/slope-models-adpkd/code/VC_12_IgAN.R", echo = FALSE)
```

After sourcing, following figure and table are saved in /data/results
folder

- Figure S4

- Table S5

Then, quit

``` r
quit()
```

## Enrichment Analysis on Whole Proteome of SC

First start with
“~/Data/slope-models-adpkd/code/SC_20_enrichmentGraph.R”. This script
will generate figures and tables.

``` r
source("~/Data/slope-models-adpkd/code/SC_20_enrichmentGraph.R", echo = FALSE)
```

After sourcing, following figures and data are saved in /data/results
folder

- Figure S2, S3

- Data S1, S2
