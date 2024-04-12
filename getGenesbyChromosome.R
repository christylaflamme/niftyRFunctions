# get a list of all the genes on X or Y
# useful code snippet
# get list of genes on a specific chromosome from bioMart
# Christy LaFlamme
# 04-12-2024

# load libraries
library(biomaRt)
library(dplyr)

# acquire mart object from human ensembl dataset
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# grab chromosome and hgnc symbol
results <- getBM(attributes = c("chromosome_name", "hgnc_symbol"), mart = mart)

# subset the results by chromosome
XY <- results %>% filter(chromosome_name == c("X", "Y"))
head(XY)

# chromosome_name hgnc_symbol
# 1               Y     TSPY24P
# 2               Y     RBMY3AP
# 3               Y     TSPY12P
# 4               Y     DUX4L19
# 5               Y      TTTY7B
# 6               Y     RBMY2HP

dim(XY) # [1] 1178    2

sessionInfo()

# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: /UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4        biomaRt_2.58.2     SeuratObject_5.0.1 sp_2.1-3          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2                   bitops_1.0-7                rlang_1.1.3                 magrittr_2.0.3              matrixStats_1.2.0           compiler_4.3.1              RSQLite_2.3.5              
# [8] DelayedMatrixStats_1.24.0   png_0.1-8                   vctrs_0.6.5                 reshape2_1.4.4              stringr_1.5.1               pkgconfig_2.0.3             crayon_1.5.2               
# [15] fastmap_1.1.1               dbplyr_2.4.0                XVector_0.42.0              utf8_1.2.4                  purrr_1.0.2                 bit_4.0.5                   zlibbioc_1.48.2            
# [22] cachem_1.0.8                beachmat_2.18.1             GenomeInfoDb_1.38.7         progress_1.2.3              blob_1.2.4                  DelayedArray_0.28.0         BiocParallel_1.36.0        
# [29] irlba_2.3.5.1               parallel_4.3.1              prettyunits_1.2.0           R6_2.5.1                    stringi_1.8.3               parallelly_1.37.1           GenomicRanges_1.54.1       
# [36] Rcpp_1.0.12                 SummarizedExperiment_1.32.0 future.apply_1.11.1         IRanges_2.36.0              Matrix_1.6-5                tidyselect_1.2.1            rstudioapi_0.15.0          
# [43] abind_1.4-5                 codetools_0.2-19            curl_5.2.1                  listenv_0.9.1               lattice_0.22-5              tibble_3.2.1                plyr_1.8.9                 
# [50] withr_3.0.0                 Biobase_2.62.0              KEGGREST_1.42.0             future_1.33.1               BiocFileCache_2.10.1        xml2_1.3.6                  Biostrings_2.70.3          
# [57] filelock_1.0.3              pillar_1.9.0                MatrixGenerics_1.14.0       stats4_4.3.1                generics_0.1.3              RCurl_1.98-1.14             S4Vectors_0.40.2           
# [64] hms_1.1.3                   ggplot2_3.5.0               sparseMatrixStats_1.14.0    munsell_0.5.0               scales_1.3.0                globals_0.16.3              glue_1.7.0                 
# [71] tools_4.3.1                 ScaledMatrix_1.10.0         locfit_1.5-9.9              dotCall64_1.1-1             XML_3.99-0.16.1             cowplot_1.1.3               grid_4.3.1                 
# [78] AnnotationDbi_1.64.1        colorspace_2.1-0            GenomeInfoDbData_1.2.11     patchwork_1.2.0             BiocSingular_1.18.0         cli_3.6.2                   rsvd_1.0.5                 
# [85] rappdirs_0.3.3              spam_2.10-0                 fansi_1.0.6                 S4Arrays_1.2.1              gtable_0.3.4                DESeq2_1.42.1               digest_0.6.35              
# [92] progressr_0.14.0            BiocGenerics_0.48.1         SparseArray_1.2.4           ggrepel_0.9.5               dqrng_0.3.2                 memoise_2.0.1               lifecycle_1.0.4            
# [99] httr_1.4.7                  PCAtools_2.14.0             bit64_4.0.5                
