---

layout: page
title: Tutorial
<!-- subtitle: An R/Bioconductor software tool to identify genes associated with latent space interactions in spatial transcriptomics. -->
hero_image: /spacemarkers/images/spacehero.jpg
<!-- hero_height: is-fullwidth -->
hero_darken: true
show_sidebar: false
hero_link: https://github.com/FertigLab/SpaceMarkers
hero_link_text: GitHub Repository

---

# Installation

You can install SpaceMarkers directly from the GitHub source.

```
install.packages(“remotes”)
remotes::install_github(“FertigLab/SpaceMarkers”, dependencies = TRUE, build_vignettes = TRUE)
```

# Importing Libraries

```
library(SpaceMarkers)
#> Warning: replacing previous import 'ape::where' by 'dplyr::where' when loading
#> 'SpaceMarkers'
#> Warning: replacing previous import 'dplyr::count' by 'matrixStats::count' when
#> loading 'SpaceMarkers'
#> Warning: replacing previous import 'rstatix::filter' by 'stats::filter' when
#> loading 'SpaceMarkers'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'SpaceMarkers'
```

# Obtaining and Formatting the Data

The data that will be used to demonstrate SpaceMarkers capabilities is a human breast cancer spatial transcriptomics dataset that comes from Visium. The CoGAPS patterns as seen in the manuscript [Atul Deshpande, Melanie Loth, et al.](https://www.biorxiv.org/content/10.1101/2022.06.02.490672v1) will also be taken from GitHub.

```
counts_url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
spatial_url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz"
cogaps_url <- "https://github.com/atuldeshpande/SpaceMarkers-paper/blob/main/CoGAPS_Analysis/BreastCancer/182967_1168993F_2_CogapsResult_5.rds?raw=true"
opt_params_path <- "optParams_breast_cancer.tsv"
```
