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

```r
install.packages(“remotes”)
remotes::install_github(“FertigLab/SpaceMarkers”, dependencies = TRUE, build_vignettes = TRUE)
```

# Importing Libraries

```r
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

```r
counts_url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
spatial_url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz"
cogaps_url <- "https://github.com/atuldeshpande/SpaceMarkers-paper/blob/main/CoGAPS_Analysis/BreastCancer/182967_1168993F_2_CogapsResult_5.rds?raw=true"
opt_params_path <- "optParams_breast_cancer.tsv"
```

# Extracting Counts Matrix

Here the counts matrix will be obtained from the h5 object in the Visium site and genes with less than 3 counts are removed from the dataset.

```r
system(paste0("wget -q ",counts_url))
counts_matrix <- load10XExpr(visiumDir = ".", h5filename = basename(counts_url))
#> Warning in asMethod(object): sparse->dense coercion: allocating vector of size
#> 1.3 GiB
good_gene_threshold <- 3
goodGenes <- rownames(counts_matrix)[apply(counts_matrix,1,function(x) sum(x>0)>=good_gene_threshold)]
counts_matrix <- counts_matrix[goodGenes,]
```

# Obtaining CoGAPS Patterns

In this example the latent features from CoGAPS will be used to identify overlapping genes with SpaceMarkers. Here the featureLoadings (cells) and samplePatterns (genes) for both the expression matrix and CoGAPS matrix need to match.

```r
system(paste0("wget -q ",cogaps_url))
cogaps_result <- readRDS(basename(cogaps_url))
features <- intersect(rownames(counts_matrix),rownames(cogaps_result@featureLoadings))
barcodes <- intersect(colnames(counts_matrix),rownames(cogaps_result@sampleFactors))
counts_matrix <- counts_matrix[features,barcodes]
cogaps_matrix <- cogaps_result@featureLoadings[features,] %*% t(cogaps_result@sampleFactors[barcodes,])
```

# Obtaining Spatial Coordinates

The spatial coordinates will also be pulled from Visium for this dataset. These are combined with the latent features to demonstrate how cells for each pattern interact in 2D space.

```r
download.file(spatial_url, basename(spatial_url))
untar(basename(spatial_url))
spCoords <- load10XCoords(visiumDir = ".")
rownames(spCoords) <- spCoords$barcode
spCoords <- spCoords[barcodes,]
spPatterns <- cbind(spCoords,cogaps_result@sampleFactors[barcodes,])
head(spPatterns)
#>                               barcode         y        x    Pattern_1
#> AAACAACGAATAGTTC-1 AAACAACGAATAGTTC-1  79.44776 163.2886 0.4676255882
#> AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1 355.49321 435.6841 0.2690758109
#> AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1  96.05857 248.8780 0.1105933860
#> AAACACCAATAACTGC-1 AAACACCAATAACTGC-1 404.91037 172.5120 0.0002508377
#> AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1 156.88474 410.5056 0.2849308848
#> AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1 316.63265 140.8859 0.1583736390
#>                       Pattern_2    Pattern_3    Pattern_4    Pattern_5
#> AAACAACGAATAGTTC-1 1.049391e-01 2.576064e-01 0.6848062277 4.747092e-02
#> AAACAAGTATCTCCCA-1 4.394425e-01 2.056469e-01 0.2921337187 1.167576e-02
#> AAACAATCTACTAGCA-1 1.148523e-02 2.309153e-01 0.4111314714 9.508318e-02
#> AAACACCAATAACTGC-1 1.685795e-01 1.223603e-01 0.0001562788 8.041928e-01
#> AAACAGAGCGACTCCT-1 1.102506e-01 9.053156e-08 0.2429406196 3.430807e-08
#> AAACAGCTTTCAGAAG-1 9.741083e-06 1.723470e-01 0.3059957027 7.167605e-01
```
