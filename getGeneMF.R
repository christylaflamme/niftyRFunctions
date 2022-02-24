###############################################################################
# getGeneMF
# 02-24-22
# Christy LaFlamme
# Purpose: takes a list of genes and adds gene functional information 
###############################################################################
# genelist: should be a list of genes (RefSeq Names)
# workDir: should be a string path to the working directory
# queryName: should be a string name for the output file
###############################################################################

# load dependencies
library(mygene)
library(data.table)
library(readxl)

getGeneMF <- function(genelist, workDir, queryName) {
  
  # query functional information
  res <- queryMany(genelist, scopes = 'symbol', fields = c('entrezgene', 'go'), species = 'human')
  res.df <- data.frame(res) 
  res.pos <- res.df[is.na(res$notfound),] # only keep rows for which annotation information was found
  
  # write a loop that accumulates MF info for all successful queries
  for (i in 1:nrow(res.pos)){ # loop over the queries
    gene <- NULL
    res.MF <- NULL
    
    gene <- res.pos$query[i] # loop over the genes
    res.MF <- res[res$query == gene, 'go.MF'][[1]] # get MF information
    res.pos$MF_summary[i] <- toString(res.MF$term) # add MF information to the dataframe
    
  }
  
  res.print <- res.pos[,c(1,9)] # print the gene name and MF information
  
  write.table(res.print, file = paste(workDir, queryName, sep = ""), quote = F, sep = "\t", row.names= F) # write to text file 
}

###############################################################################
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.4
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] readxl_1.3.1           data.table_1.14.2      mygene_1.28.0          GenomicFeatures_1.44.2 AnnotationDbi_1.54.1   Biobase_2.52.0         GenomicRanges_1.44.0  
# [8] GenomeInfoDb_1.28.4    IRanges_2.26.0         S4Vectors_0.30.2       BiocGenerics_0.38.0   
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7                matrixStats_0.61.0          bit64_4.0.5                 filelock_1.0.2              RColorBrewer_1.1-2          progress_1.2.2             
# [7] httr_1.4.2                  tools_4.1.1                 backports_1.3.0             utf8_1.2.2                  R6_2.5.1                    rpart_4.1-15               
# [13] Hmisc_4.6-0                 DBI_1.1.1                   colorspace_2.0-2            nnet_7.3-16                 tidyselect_1.1.1            gridExtra_2.3              
# [19] prettyunits_1.1.1           chron_2.3-56                bit_4.0.4                   curl_4.3.2                  compiler_4.1.1              htmlTable_2.3.0            
# [25] xml2_1.3.2                  DelayedArray_0.18.0         rtracklayer_1.52.1          scales_1.1.1                checkmate_2.0.0             rappdirs_0.3.3             
# [31] stringr_1.4.0               digest_0.6.28               Rsamtools_2.8.0             foreign_0.8-81              XVector_0.32.0              htmltools_0.5.2            
# [37] base64enc_0.1-3             jpeg_0.1-9                  pkgconfig_2.0.3             MatrixGenerics_1.4.3        dbplyr_2.1.1                fastmap_1.1.0              
# [43] htmlwidgets_1.5.4           rlang_0.4.12                rstudioapi_0.13             RSQLite_2.2.8               BiocIO_1.2.0                generics_0.1.1             
# [49] jsonlite_1.7.2              BiocParallel_1.26.2         dplyr_1.0.7                 RCurl_1.98-1.5              magrittr_2.0.1              GenomeInfoDbData_1.2.6     
# [55] Formula_1.2-4               Matrix_1.3-4                Rcpp_1.0.7                  munsell_0.5.0               fansi_0.5.0                 proto_1.0.0                
# [61] sqldf_0.4-11                lifecycle_1.0.1             stringi_1.7.5               yaml_2.2.1                  SummarizedExperiment_1.22.0 zlibbioc_1.38.0            
# [67] plyr_1.8.6                  BiocFileCache_2.0.0         grid_4.1.1                  blob_1.2.2                  crayon_1.4.2                lattice_0.20-45            
# [73] Biostrings_2.60.2           splines_4.1.1               hms_1.1.1                   KEGGREST_1.32.0             knitr_1.36                  pillar_1.6.4               
# [79] rjson_0.2.20                biomaRt_2.48.3              XML_3.99-0.8                glue_1.5.0                  latticeExtra_0.6-29         png_0.1-7                  
# [85] vctrs_0.3.8                 cellranger_1.1.0            gtable_0.3.0                purrr_0.3.4                 gsubfn_0.7                  assertthat_0.2.1           
# [91] cachem_1.0.6                ggplot2_3.3.5               xfun_0.28                   restfulr_0.0.13             survival_3.2-13             tibble_3.1.6               
# [97] GenomicAlignments_1.28.0    memoise_2.0.0               cluster_2.1.2               ellipsis_0.3.2   








