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
