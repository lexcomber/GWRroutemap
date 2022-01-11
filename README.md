# A route map for successful applications of Geographically Weighted Regression

Alexis Comber<sup>1*</sup>, Chris Brunsdon<sup>2</sup>, Martin Charlton<sup>2</sup>, Guanpeng Dong<sup>3</sup>, Rich Harris<sup>4</sup>, Binbin Lu<sup>5</sup>, Yihe Lü<sup>6</sup>, Daisuke Murakami<sup>7</sup>, Tomoki Nakaya<sup>8</sup>, Yunqiang Wang<sup>9</sup>, Paul Harris<sup>10</sup>

<sup>1</sup> School of Geography, University of Leeds, Leeds, UK.\
<sup>2</sup> National Centre for Geocomputation, Maynooth University, Maynooth, Ireland.\
<sup>3</sup> Key Research Institute of Yellow River Civilization and Sustainable Development, Henan University, Kaifeng, China.\
<sup>4</sup> School of Geographical Sciences, University of Bristol, Bristol, UK.\
<sup>5</sup> School of Remote Sensing and Information Engineering, Wuhan University, Wuhan, China.\
<sup>6</sup> State Key Laboratory of Urban and Regional Ecology, Research Center for Eco-Environmental Sciences, Chinese Academy of Sciences; Joint Center for Global Change Studies; University of Chinese Academy of Sciences, Beijing, China.\
<sup>7</sup> Department of Statistical Data Science, Institute of Statistical Mathematics, Tachikawa, Japan.\
<sup>8</sup> Graduate School of Environmental Studies, Tohoku University, Sendai, Japan.\
<sup>9</sup> State Key Laboratory of Loess and Quaternary Geology, Institute of Earth Environment, Chinese Academy of Sciences, Xi’an, China.\
<sup>10</sup> Sustainable Agriculture Sciences North Wyke, Rothamsted Research, Okehampton, UK.

<sup>*</sup> contact author: a.comber@leeds.ac.uk

## Abstract

Geographically Weighted Regression (GWR) is increasingly used in spatial analyses of social and environmental data. It allows spatial heterogeneities in processes and relationships to be investigated through a series of local regression models rather than a single global one. Standard GWR assumes that relationships between the response and predictor variables operate at the same spatial scale, which is frequently not the case. To address this, several GWR variants have been proposed. This paper describes a route map to decide whether to use a GWR model or not, and if so which of three core variants to apply: a standard GWR, a mixed GWR or a multiscale GWR (MS-GWR). The route map comprises 3 primary steps that should always be undertaken: (1) a basic linear regression, (2) a MS-GWR, and (3) investigations of the results of these in order to decide whether to use a GWR approach, and if so for determining the appropriate GWR variant. The paper also highlights the importance of investigating a number of secondary issues at global and local scales including collinearity, the influence of outliers, and dependent error terms. Code and data for the case study used to illustrate the route map are provided.

**Keywords**: Spatially varying coefficient model; non-stationarity; spatial heterogeneity; autocorrelation; regression

The paper is published in Geographical Analysis (Open Access at https://onlinelibrary.wiley.com/doi/epdf/10.1111/gean.12316) and the draws heavily from an earier version now archived at https://arxiv.org/abs/2004.06070. 


## Code
To run the analysis in this paper you should download the the R script `GWR_route_map_git_2022.R` and install the packages. Package and other info is below. The data are pulled from the GitHub site. The code recreates the results as the same sequence in the paper. 

If you have any problems with data / code / versions etc please contact Lex Comber at the email above.

> sessionInfo()
`R version 4.0.4 (2021-02-15)`
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] deldir_1.0-6        nlme_3.1-152        gridExtra_2.3       broom_0.7.6         car_3.0-10          carData_3.0-4      
 [7] OpenStreetMap_0.3.4 repmis_0.5          forcats_0.5.1       stringr_1.4.0       dplyr_1.0.5         purrr_0.3.4        
[13] readr_1.4.0         tidyr_1.1.3         tibble_3.1.1        ggplot2_3.3.5       tidyverse_1.3.1     GWmodel_2.2-8      
[19] spatialreg_1.1-5    Matrix_1.3-2        Rcpp_1.0.7          robustbase_0.93-7   spdep_1.1-7         sf_0.9-8           
[25] spData_2.0.1        rgdal_1.5-23        raster_3.5-2        GISTools_0.7-4      rgeos_0.5-5         MASS_7.3-53.1      
[31] RColorBrewer_1.1-2  maptools_1.1-2      sp_1.4-6           

loaded via a namespace (and not attached):
 [1] colorspace_2.0-0   ellipsis_0.3.1     class_7.3-18       rio_0.5.26         fs_1.5.0           rstudioapi_0.13   
 [7] proxy_0.4-26       farver_2.1.0       fansi_0.4.2        lubridate_1.7.10   xml2_1.3.2         codetools_0.2-18  
[13] splines_4.0.4      R.methodsS3_1.8.1  knitr_1.32         jsonlite_1.7.2     rJava_1.0-5        dbplyr_2.1.1      
[19] R.oo_1.24.0        compiler_4.0.4     httr_1.4.2         backports_1.2.1    assertthat_0.2.1   cli_2.4.0         
[25] tools_4.0.4        coda_0.19-4        gtable_0.3.0       glue_1.4.2         gmodels_2.18.1     cellranger_1.1.0  
[31] vctrs_0.3.7        gdata_2.18.0       xfun_0.22          ps_1.6.0           openxlsx_4.2.3     rvest_1.0.0       
[37] lifecycle_1.0.0    gtools_3.9.2       terra_1.4-11       DEoptimR_1.0-9     LearnBayes_2.15.1  zoo_1.8-9         
[43] scales_1.1.1       hms_1.0.0          parallel_4.0.4     expm_0.999-6       curl_4.3           stringi_1.5.3     
[49] e1071_1.7-9        boot_1.3-27        zip_2.1.1          intervals_0.15.2   rlang_0.4.10       pkgconfig_2.0.3   
[55] lattice_0.20-41    labeling_0.4.2     tidyselect_1.1.0   plyr_1.8.6         magrittr_2.0.1     R6_2.5.0          
[61] generics_0.1.0     DBI_1.1.1          pillar_1.6.0       haven_2.4.0        foreign_0.8-81     withr_2.4.2       
[67] units_0.7-1        xts_0.12.1         abind_1.4-5        spacetime_1.2-5    modelr_0.1.8       crayon_1.4.1      
[73] KernSmooth_2.23-18 utf8_1.2.1         grid_4.0.4         readxl_1.3.1       data.table_1.14.0  FNN_1.1.3         
[79] reprex_2.0.0       digest_0.6.27      classInt_0.4-3     R.cache_0.14.0     R.utils_2.10.1     munsell_0.5.0    
