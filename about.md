---

layout: page
title: About
<!-- subtitle: An R/Bioconductor software tool to identify genes associated with latent space interactions in spatial transcriptomics. -->
hero_image: /SpaceMarkersGuide/images/spacehero.jpg
<!-- hero_height: is-fullwidth -->
hero_darken: true
show_sidebar: false
<!-- hero_link: https://github.com/FertigLab/SpaceMarkers
hero_link_text: GitHub Repository -->

---

<img src="/SpaceMarkersGuide/images/SpaceMarkersHexWhite250.png" align="left" style="margin: 0px 20px 0px 0px;" />**SpaceMarkers** is an R package which identifies genes related to cell-cell interactions in **spatial transcriptomics (ST) data**. SpaceMarkers uses latent space analysis to identify spatially overlapping latent features associated with different cellular signatures, inferring the genes associated with their interaction. This method relies on estimating spatially resolved linear latent features that represent cellular signatures in the ST data. These features are characterized by continuous weights corresponding to each spatial coordinate in the ST data, and the patterns in the ST data are denoted by these continuous weights. The **input** to the SpaceMarkers algorithm is the ST data matrix and spatially resolved patterns learned through latent space analysis, and the **output** is a list of genes associated with the interaction between each pair of spatially overlapping patterns.

The **first stage** of SpaceMarkers involves identifying each pattern’s region of influence and subsequently the region of pattern interaction. When two patterns have overlapping influence in the same region of the tissue, we assume an interaction between these patterns in this interaction region. The **second stage** ranks genes exhibiting higher activity levels in the interaction region relative to regions with exclusive influence from each pattern.

SpaceMarkers provides a “**differential expression**” (DE) mode to quantify genes with enhanced expression in a region with overlapping influence from two patterns compared to regions with exclusive influence from individual patterns. It also provides a “**residual**” mode, which identifies genes that have significantly higher residual error between the original ST data and its estimated fit from the CoGAPS model in the region with overlapping influence from two patterns when compared to the regions with exclusive influence from each pattern. This approach enables the detection of effects of intercellular interaction after taking into account variations in cell populations. In both modes, the algorithm uses a **non-parametric Kruskal-Wallis statistical test** in the three regions of interest followed by **posthoc analysis** to identify these genes, which make up the SpaceMarkers output.

While the examples shown here use CoGAPS, SpaceMarkers is generally applicable to latent features from any of a number of linear latent feature factorization approaches available in the literature.
