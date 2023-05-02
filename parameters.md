---

layout: page
title: Parameters
<!-- subtitle: An R/Bioconductor software tool to identify genes associated with latent space interactions in spatial transcriptomics. -->
hero_image: /spacemarkers/images/spacehero.jpg
<!-- hero_height: is-fullwidth -->
hero_darken: true
show_sidebar: false
<!-- hero_link: https://github.com/FertigLab/SpaceMarkers
hero_link_text: GitHub Repository -->

---

# load10XExpr() Arguments

<table>
<thead>
<tr class="header">
<th align="left">Argument</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">visiumDir</td>
<td align="left">A string path to the h5 file with expression
information</td>
</tr>
<tr class="even">
<td align="left">h5filename</td>
<td align="left">A string of the name of the h5 file in the
directory</td>
</tr>
</tbody>
</table>

# load10XCoords() Arguments

<table>
<colgroup>
<col width="4%" />
<col width="95%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Argument</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">visiumDir</td>
<td align="left">A string path to the location of the folder containing
the spatial coordinates. The folder in your visiumDir must be named
‘spatial’ and must contain files ‘scalefactors_json.json’ and
‘tissue_positions_list.csv’.</td>
</tr>
<tr class="even">
<td align="left">resolution</td>
<td align="left">A string specifying which values to look for in the
.json object. Can be either lowres or highres.</td>
</tr>
</tbody>
</table>

# getSpatialParameters() Arguments

<table>
<colgroup>
<col width="8%" />
<col width="91%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Argument</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">spatialPatterns</td>
<td align="left">A data frame that contains the spatial coordinates for
each cell type. The column names must include ‘x’ and ‘y’ as well as a
set of numbered columns named ‘Pattern_1…..N’.</td>
</tr>
</tbody>
</table>

# getInteractingGenes() Arguments

<table>
<colgroup>
<col width="8%" />
<col width="91%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Argument</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">data</td>
<td align="left">A data frame of expression information with rows being
the features/genes and columns being the samples/cells.</td>
</tr>
<tr class="even">
<td align="left">reconstruction</td>
<td align="left">A data frame of features (rows) and samples (columns)
constructed from the information of a latent feature method such as
CoGAPs or STdeconvolve. NULL if ‘DE’ mode is specified</td>
</tr>
<tr class="odd">
<td align="left">optParams</td>
<td align="left">A data frame that for each pattern has the sigmaOpts -
the optimal width of the gaussian distribution and the thresOpt -
outlier threshold around the set of spots.</td>
</tr>
<tr class="even">
<td align="left">spatialPatterns</td>
<td align="left">A data frame that contains the spatial coordinates for
each cell type. The column names must include ‘x’ and ‘y’ as well as a
set of numbered columns named ‘Pattern_1…..N’.</td>
</tr>
<tr class="odd">
<td align="left">refPattern</td>
<td align="left">A string of the pattern you want to use to compare to
the other patterns in the latent feature space</td>
</tr>
<tr class="even">
<td align="left">mode</td>
<td align="left">A string specifying either ‘residual’ mode or ‘DE’ mode
for finding interacting genes</td>
</tr>
<tr class="odd">
<td align="left">minOverlap</td>
<td align="left">a number that specifies the minimum overlap between
genes in two patterns to be considered for the statistical tests. The
default is 50.</td>
</tr>
<tr class="even">
<td align="left">hotspotRegions</td>
<td align="left">a vector that specifies the patterns to compare to the
‘refPattern’. The default is NULL which indicates that all patterns
would be compared to the ‘refPattern’.</td>
</tr>
</tbody>
</table>
