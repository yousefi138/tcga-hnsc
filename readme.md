# Using genomic data to predict Head and Neck cancer outcomes

Background and instructions presentation: [pdf](slides.pdf) [pptx](slides.pptx)

## System requirements

To successfully run the pipeline below and download the TCGA data
you'll need to have access to a bash/unix command line with the R statistical
programming language installed. 

No R programming should be required, as all tasks are coordinated by
bash shell scripts.

However, a few R packages must be installed prior to implementing the 
pipeline. To install these, run the below command:

```
bash install-packages.sh
```

## Downloading and preparing the dataset

Considered using the TCGAbiolinks R package to download TCGA data
but ran into installation problems.
Instead, I compiled a list of files to download from the GDAC website:
http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128
See `files.csv` in this directory.

See [variables.txt](variables.txt)) for full list of directories of where
data are downloaded, cleaned, and stored. 

All data generated are within the `PROJECT_DIR` which I've set to be `hnsc-project`
but which can be freely changed to preference. 

Data files will be downloaded to the `FILES_DIR` directory by using the 
below command:

```
bash download-data.sh files.csv
```

Clinical outcome data has been cleaned up as part of the
PanCancer Atlas project
(https://gdc.cancer.gov/about-data/publications/pancanatlas).

> Liu J, Lichtenberg T, Hoadley KA, et al. An Integrated TCGA Pan-Cancer
> Clinical Data Resource to Drive High-Quality Survival Outcome
> Analytics. Cell. 2018;173(2):400-416.e11. doi:10.1016/j.cell.2018.02.052

This publication cautions against using overall survival as an outcome
because the follow-up isn't long enough.
Recommends progression-free interval (PFI) or
disease-free interval (DFI).
PFI and DFI are available in Supplementary Table 1
(https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81).
The table is downloaded to the `PAN_CANCER_DIR` directory
using the following script.

```
bash download-pan-cancer-clinical.sh
```

The dataset will be generated
from the downloaded files to the `FULL_DIR` directory
using the followin script.

```
bash extract-data.sh
```

This dataset is then split into protein and methylation subsets
in the `PROTEIN_DIR` and `METHYLATION_DIR` directories, respectively.

```
bash split-data.sh
```

Prepare annotations of the data in order to link between
the different types. Creates three annotation files:
- `METHYLATION_DIR/methylation-annotation.txt` Links CpG sites to genes.
- `PROTEIN_DIR/protein-annotation.txt` Links proteins to genes.

```
bash annotate-data.sh
```



