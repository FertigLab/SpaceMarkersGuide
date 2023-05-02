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

## Extracting Counts Matrix

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

## Obtaining CoGAPS Patterns

In this example the latent features from CoGAPS will be used to identify overlapping genes with SpaceMarkers. Here the featureLoadings (cells) and samplePatterns (genes) for both the expression matrix and CoGAPS matrix need to match.

```r
system(paste0("wget -q ",cogaps_url))
cogaps_result <- readRDS(basename(cogaps_url))
features <- intersect(rownames(counts_matrix),rownames(cogaps_result@featureLoadings))
barcodes <- intersect(colnames(counts_matrix),rownames(cogaps_result@sampleFactors))
counts_matrix <- counts_matrix[features,barcodes]
cogaps_matrix <- cogaps_result@featureLoadings[features,] %*% t(cogaps_result@sampleFactors[barcodes,])
```

## Obtaining Spatial Coordinates

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

# Executing SpaceMarkers

## SpaceMarker Modes

SpaceMarkers can operate in ‘residual’ mode or ‘DE’ (Differential Expression mode). In an ideal world the overlapping patterns identified by SpaceMarkers would be a homogeneous population of cells and the relationship between them would be linear. However, due to confounding effects of variations in cell density and common cell types in any given region, this is not always true.

To account for these confounding effects, the ‘residual’ mode compares the feature interactions between the expression matrix and the reconstructed latent space matrix. The features with the highest residual error are reported. The genes are then classified according to regions of overlapping vs exclusive influence. The default mode is ‘residual’ mode.

Suppose the feature (gene) information is not readily available and only the sample (cells) latent feature patterns with P-values are available? This is the advantage of ‘DE’ mode. Where residual mode assesses the non-linear effects that may arise from confounding variables, ‘DE’ mode assesses simple linear interactions between patterns. DE mode also compares genes from regions of overlapping vs exclusive influence but does not consider residuals from the expression matrix as there is no matrix reconstruction with the latent feature matrix.

To demonstrate SpaceMarkers we will be looking at Pattern_1 from CoGAPS. Pattern_1 was identified as mainly an immune driven pattern. Pattern_1 is the default setting if no Pattern_n preference is given by the user.

```r
SpaceMarkersMode = "residual" 
SpaceMarkersRefPattern = "Pattern_1" 
```

## Residual Mode

SpaceMarkers identifies regions of influence using a gaussian kernel outlier based model. The reference pattern (Pattern_1 in this case) is used as the prior for this model. SpaceMarkers then identifies where the regions of influence are interacting from each of the other patterns as well as where they are mutually exclusive.

getSpatialParameters: This function identifies the optimal width of the gaussian distribution (sigmaOpt) as well as the outlier threshold around the set of spots (thresOpt) for each pattern.These parameters minimize the spatial autocorrelation residuals of the spots in the regions of influence. This function can take a while so the file would be read in from the data folder for this tutorial.

```r
#Takes approximately 12 minutes
#optParams <- getSpatialParameters(spPatterns)
optParams <- as.matrix(read.csv(opt_params_path, header = TRUE, sep = "\t"))
head(optParams)
#>           Pattern_1 Pattern_2 Pattern_3 Pattern_4 Pattern_5
#> sigmaOpt        7.0       7.8       6.0       7.2       6.4
#> threshOpt       2.1       1.3       2.6       1.5       1.2
```

getInteractingGenes: This function identifies the regions of influence and interaction as well as the genes associated with these regions. A non-parametric Kruskal-Wallis test is used to identify genes statistically significant genes in any one region of influence and a post hoc Dunn’s Test is used for analysis of genes between regions. If ‘residual’ mode is selected the user must provide a reconstructed matrix from the latent feature matrix. This is passed to the ‘reconstruction’ argument and can be left as NULL for ‘DE’ mode. The ‘data’ parameter is the original expression matrix. The ‘spatialPatterns’ argument takes a matrix with the spatial coordinates of each cell as well as the patterns. The spatial coordinate columns must have the labels x and y.

The output is a list of data frames with information about the interacting genes of the refPattern and each pattern from the CoGAPS matrix (interacting_genes object). There is also a data frame with all of the regions of influence for any two of patterns (the hotspotRegions object).

For the ‘interacting_genes’ data frames, the first column is the list of genes and the second column says whether the genes are statistically significant in the refPattern only, the interacting pattern only, or both. The remaining columns are statistics for the Kruskal-Wallis test and the post hoc Dunn’s test.

```r
print(head(SpaceMarkers$interacting_genes[[1]]))
#>                  Gene Pattern_1 x Pattern_2 KW.obs.tot KW.obs.groups KW.df
#> CD52             CD52                vsBoth       2061             3     2
#> CORO1A         CORO1A                vsBoth       2061             3     2
#> AC239799.2 AC239799.2                vsBoth       2061             3     2
#> AC015987.1 AC015987.1                vsBoth       2061             3     2
#> IGLC3           IGLC3                vsBoth       2061             3     2
#> HLA-DQB1     HLA-DQB1                vsBoth       2061             3     2
#>            KW.statistic    KW.pvalue     KW.p.adj Dunn.zP1_Int Dunn.zP2_Int
#> CD52           136.2669 2.570528e-30 1.026239e-29    -2.430175    -9.878900
#> CORO1A         117.7495 2.697858e-26 9.927205e-26    -2.470913    -9.296994
#> AC239799.2     113.3366 2.450616e-25 8.839040e-25    -9.423240    -3.019372
#> AC015987.1     109.6398 1.556057e-24 5.498136e-24    -2.532100    -9.048589
#> IGLC3          106.8834 6.174019e-24 2.154040e-23    -5.576924   -10.128681
#> HLA-DQB1       105.1078 1.500154e-23 5.192243e-23    -2.371747    -8.803429
#>            Dunn.zP2_P1 Dunn.pval_1_Int Dunn.pval_2_Int Dunn.pval_2_1
#> CD52         -9.036216    1.509155e-02    5.139376e-23  1.621886e-19
#> CORO1A       -8.260723    1.347686e-02    1.444727e-20  1.447929e-16
#> AC239799.2    8.719814    4.373794e-21    2.532992e-03  2.786595e-18
#> AC015987.1   -7.871751    1.133817e-02    1.448270e-19  3.497118e-15
#> IGLC3        -5.183776    2.448092e-08    4.121719e-24  2.174382e-07
#> HLA-DQB1     -7.779831    1.770422e-02    1.326984e-18  7.262134e-15
#>            Dunn.pval_1_Int.adj Dunn.pval_2_Int.adj Dunn.pval_2_1.adj
#> CD52              4.313251e-02        5.420730e-22      1.558207e-18
#> CORO1A            3.883638e-02        1.430031e-19      1.247539e-15
#> AC239799.2        4.393647e-20        8.375101e-03      2.568639e-17
#> AC015987.1        3.315488e-02        1.394354e-18      2.846410e-14
#> IGLC3             1.396007e-07        4.463128e-23      1.157326e-06
#> HLA-DQB1          4.993202e-02        1.237561e-17      5.832096e-14
print(head(SpaceMarkers$hotspotRegions))
#>      Pattern_1   Pattern_2   Pattern_3 Pattern_4   Pattern_5  
#> [1,] "Pattern_1" NA          NA        "Pattern_4" NA         
#> [2,] NA          "Pattern_2" NA        "Pattern_4" NA         
#> [3,] NA          NA          NA        "Pattern_4" NA         
#> [4,] NA          NA          NA        NA          "Pattern_5"
#> [5,] "Pattern_1" "Pattern_2" NA        NA          NA         
#> [6,] NA          NA          NA        NA          "Pattern_5"
```

## DE Mode

As described previously ‘DE’ mode only requires the counts matrix and spatial patterns and not the reconstructed CoGAPS matrix. It identifies simpler molecular interactions between regions.

```r

SpaceMarkersMode = "DE"
SpaceMarkers_DE <- getInteractingGenes(data = counts_matrix, reconstruction = NULL, optParams = optParams,spatialPatterns = spPatterns, refPattern = SpaceMarkersRefPattern, mode = SpaceMarkersMode)
#> [1] "Using user provided optParams."
#> Warning in matrixTests::row_kruskalwallis(x = as.matrix(testMat), g = region):
#> 2831 columns dropped due to missing group information
#> Warning: matrixTests::row_kruskalwallis: 209 of the rows were essentially constant.
#> First occurrence at row 3
#> Warning in matrixTests::row_kruskalwallis(x = as.matrix(testMat), g = region):
#> 2979 columns dropped due to missing group information
#> Warning: matrixTests::row_kruskalwallis: 128 of the rows were essentially constant.
#> First occurrence at row 67
#> Warning in matrixTests::row_kruskalwallis(x = as.matrix(testMat), g = region):
#> 2671 columns dropped due to missing group information
#> Warning: matrixTests::row_kruskalwallis: 131 of the rows were essentially constant.
#> First occurrence at row 3
#> Warning in matrixTests::row_kruskalwallis(x = as.matrix(testMat), g = region):
#> 2502 columns dropped due to missing group information
#> Warning: matrixTests::row_kruskalwallis: 186 of the rows were essentially constant.
#> First occurrence at row 176
```

The warnings above are due to the nature of the ‘sparse’ data being used. Comparing two cells from the two patterns with identical information is redundant as SpaceMarkers is identifying statistically different expression for interactions exclusive to either of the two patterns and a region that is due to interaction between the given two patterns. Also, if there are too many zeros in the genes (rows) of those regions, the columns are dropped as there is nothing to compare in the kruskal wallis test.

## Differences between Residual Mode and DE Mode

To highlight the differences between residual mode and DE mode, the interaction between Pattern1 (immune) and Pattern5 (invasive) will be assessed.

One of the first things to notice is the difference in the number of genes identified between the two modes.

```r
residual_p1_p2 <- SpaceMarkers$interacting_genes[[4]]
DE_p1_p2 <- SpaceMarkers_DE$interacting_genes[[4]]
```
```r
#> [1] "Residual mode identified 445 interacting genes, while DE mode identified 2715 interacting genes"
```

DE mode identified more genes because it does not consider noise that may be associated with common cell types in the interacting regions while residual mode considers these and other confounding variables by taking the residuals between the counts and reconstructed latent feature matrix.

The next analysis will show where the top genes rank in each mode’s list if they are identified at all. A function was created that will take the top 20 genes of a reference list of genes and compares it to the entire list of a second list of genes. The return object is a data frame of the gene, the name of each list and the ranking of each gene as compared to the reference list. If there is no gene identified in the second list compared to the reference it is classified as NA.

```r
res_to_DE <- compare_genes( head(residual_p1_p2$Gene, n = 20) ,DE_p1_p2$Gene,ref_name = "residual",list2_name = "DE" )

DE_to_res <- compare_genes(head(DE_p1_p2$Gene, n = 20),residual_p1_p2$Gene,ref_name = "DE",list2_name = "residual" )
```
```r
res_to_DE
#>          Gene residual_Rank DE_Rank
#> 7      CXCL14             1      NA
#> 2       APOC1             2      17
#> 12       IGHE             3       6
#> 20   Z93241.1             4      NA
#> 14       PMM2             5      NA
#> 1     ADORA2A             6      NA
#> 4        C1QB             7      NA
#> 17 SUCLG2-AS1             8      NA
#> 15     SIRLNT             9      NA
#> 8       GDF15            10      NA
#> 5        C1QC            11      21
#> 16    SLC39A4            12      NA
#> 3      B3GAT2            13      NA
#> 18       TGM2            14      24
#> 11      IFI30            15     104
#> 19     TUBA1C            16      NA
#> 6     CREB3L4            17      NA
#> 9    HLA-DRB5            18      NA
#> 10    HNRNPAB            19       9
#> 13       MAL2            20      NA
```

Here we identify the top 20 genes in ‘residual’ mode and their corresponding ranking in DE mode. IGHE, APOC1,C1QC, TGM2 and HNRNPAB are ranked high in both ‘DE’ and ‘residual’ mode. While 14 genes that are ranked high in ‘residual’ mode are not identified at all in ‘DE’ mode

```r
DE_to_res
#>       Gene DE_Rank residual_Rank
#> 11   NUPR1       1            NA
#> 19  TUBA1B       2            NA
#> 2     CD74       3            NA
#> 13   PSMB4       4           121
#> 4     CYC1       5            66
#> 7     IGHE       6             3
#> 15    RPSA       7            NA
#> 3      CLU       8            25
#> 6  HNRNPAB       9            19
#> 18  SQSTM1      10            NA
#> 16 SELENOM      11            NA
#> 20  TYROBP      12            24
#> 14   RPLP1      13            NA
#> 8    MGST3      14           151
#> 12  POLR2L      15            NA
#> 9    NINJ1      16            NA
#> 1    APOC1      17             2
#> 5      FN1      18            NA
#> 17   SERF2      19            NA
#> 10    NQO1      20            96
```

In addition to IGHE, HNRNPAB, and APOC1, CLU and TYROBP are also ranked high in both methods. There were 11 genes that were interacting in ‘DE’ mode but not ‘residual’ mode.

There is some agreement with interacting genes between the two methods but there are also quite a few differences. Therefore, the selected mode can significantly impact the downstream results and should be taken into consideration based on the specific biological question being answered and the data available.

<strong>[SpaceMarkers Parameters](/parameters/)</strong>

```r
sessionInfo()
#> R version 4.2.3 (2023-03-15)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Catalina 10.15.7
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] SpaceMarkers_0.1.0     spatstat_3.0-3         spatstat.linnet_3.0-4 
#>  [4] spatstat.model_3.1-2   rpart_4.1.19           spatstat.explore_3.0-6
#>  [7] nlme_3.1-162           spatstat.random_3.1-3  spatstat.geom_3.0-6   
#> [10] spatstat.data_3.0-0    rstatix_0.7.1          pracma_2.4.2          
#> [13] matrixTests_0.1.9.1    matrixStats_0.63.0     jsonlite_1.8.4        
#> [16] hdf5r_1.3.8            dplyr_1.1.0            ape_5.6-2             
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.10           lattice_0.20-45       tidyr_1.3.0          
#>  [4] deldir_1.0-6          digest_0.6.31         utf8_1.2.2           
#>  [7] spatstat.core_2.4-4   plyr_1.8.8            R6_2.5.1             
#> [10] backports_1.4.1       evaluate_0.20         ggplot2_3.4.0        
#> [13] tensor_1.5            pillar_1.8.1          rlang_1.0.6          
#> [16] rstudioapi_0.14       car_3.1-1             jquerylib_0.1.4      
#> [19] Matrix_1.5-3          goftest_1.2-3         qvalue_2.28.0        
#> [22] rmarkdown_2.20        splines_4.2.3         stringr_1.5.0        
#> [25] munsell_0.5.0         polyclip_1.10-4       bit_4.0.5            
#> [28] broom_1.0.3           compiler_4.2.3        xfun_0.36            
#> [31] pkgconfig_2.0.3       mgcv_1.8-42           htmltools_0.5.4      
#> [34] tidyselect_1.2.0      tibble_3.1.8          fansi_1.0.4          
#> [37] grid_4.2.3            gtable_0.3.1          lifecycle_1.0.3      
#> [40] magrittr_2.0.3        scales_1.2.1          stringi_1.7.12       
#> [43] cli_3.6.0             cachem_1.0.6          carData_3.0-5        
#> [46] reshape2_1.4.4        bslib_0.4.2           spatstat.utils_3.0-1 
#> [49] generics_0.1.3        vctrs_0.5.2           tools_4.2.3          
#> [52] bit64_4.0.5           glue_1.6.2            purrr_1.0.1          
#> [55] abind_1.4-5           parallel_4.2.3        fastmap_1.1.0        
#> [58] yaml_2.3.7            colorspace_2.1-0      spatstat.sparse_3.0-0
#> [61] knitr_1.42            sass_0.4.5
```
