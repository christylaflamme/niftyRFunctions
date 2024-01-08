######################A FUNCTION TO ASSESS EPIGENETIC AGE FOR MULTIPLE SAMPLES######################
######################
# required R packages: Sesame, data.table
library(sesame)
library(data.table)
######################
# @param: betas: a dataframe of beta values for the samples to predict epigenetic age
#         row.names should be equal to the probeIDs and the number of columns in the dataframe should be equal the exact number of samples
# @param: phenoData: a dataframe of meta data for the samples
#         "Sample_ID" column should contain the full sentrix ID and row/column numbers (e.g. "205707890075_R06C01")
#         the number of rows should be equal to the exact number of samples in the betas dataframe
######################
######################
# for sesame <predictAge> function to work, need to download and load models (see https://github.com/zhou-lab/InfiniumAnnotationV1)
# documentation: https://zhou-lab.github.io/sesame/v1.18/inferences.html#Age__Epigenetic_Clock

getAges <- function(betas, phenoData, model) {
  
  if (!ncol(betas) == nrow(phenoData)) {
    print("Number of columns in betas does not equal number of rows in phenoData. Please re-check the input data.")
  }
  
  else {
    
    # start a running column to input the epigenetic age into the phenoData
    phenoData$Epigenetic_Age <- "In_progress"
    
    for (i in 1:ncol(betas)) {
      
      # set up objects for each sample
      n.betas <- NA
      n.betas <- betas[,i] # loop over each sample of the beta value matrix
      sampleID <- toString(colnames(betas)[i]) # store the sampleID information 
      names(n.betas) <- row.names(betas)
      print(paste("Betas are now loaded for", sampleID, sep = " "))
      
      # perform epigenetic age calculation
      age <- NA
      age <- predictAge(n.betas, model) 
      phenoData[phenoData$Sample_ID == sampleID,]$Epigenetic_Age <- age # add the age information to the phenoData object where the sample matches the Sample_ID
      print(paste("The age calculation for", sampleID, "is", age, sep = " "))
      
    }
    
  }
  
  return(phenoData)
}

# example running the function
# load the betas
betas <- as.data.frame(fread("/path/to/betas.txt"))
row.names(betas) <- betas[,1] # make the probeID the row.names if first column
betas <- betas[,c(5:46)] # remove the chromosomal information if stored in betas file
phenoData <- readRDS("/path/to/phenoData.txt")

# import the model
model <- readRDS("/path/to/Clock_Horvath353.rds")

# run the function
phenoData.getAges <- getAges(betas, phenoData, model)

# save the data
write.table(phenoData.getAges, file = "/path/to/output/file.txt", sep = "\t", quote = F, row.names = FALSE)
######################
# model files for EPICv2 probeIDs corresponding to "Anno/HM450/Clock_Horvath353.rds" and "Anno/HM450/Clock_SkinBlood.rds" may be found in <zhou-lab_epigenetic_age_models>
######################
######################
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: /UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] data.table_1.14.10          sesame_1.20.0               sesameData_1.20.0           ExperimentHub_2.10.0        AnnotationHub_3.10.0        BiocFileCache_2.10.1        dbplyr_2.4.0               
# [8] minfi_1.48.0                bumphunter_1.44.0           locfit_1.5-9.8              iterators_1.0.14            foreach_1.5.2               Biostrings_2.70.1           XVector_0.42.0             
# [15] SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0       matrixStats_1.2.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.5         IRanges_2.36.0             
# [22] S4Vectors_0.40.2            BiocGenerics_0.48.1        
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3            rstudioapi_0.15.0             magrittr_2.0.3                GenomicFeatures_1.54.1        BiocIO_1.12.0                 zlibbioc_1.48.0               vctrs_0.6.5                  
# [8] multtest_2.58.0               memoise_2.0.1                 Rsamtools_2.18.0              DelayedMatrixStats_1.24.0     RCurl_1.98-1.13               askpass_1.2.0                 htmltools_0.5.7              
# [15] S4Arrays_1.2.0                progress_1.2.3                curl_5.2.0                    Rhdf5lib_1.24.1               SparseArray_1.2.3             rhdf5_2.46.1                  nor1mix_1.3-2                
# [22] plyr_1.8.9                    cachem_1.0.8                  GenomicAlignments_1.38.0      mime_0.12                     lifecycle_1.0.4               pkgconfig_2.0.3               Matrix_1.6-4                 
# [29] R6_2.5.1                      fastmap_1.1.1                 GenomeInfoDbData_1.2.11       shiny_1.8.0                   digest_0.6.33                 colorspace_2.1-0              siggenes_1.76.0              
# [36] reshape_0.8.9                 AnnotationDbi_1.64.1          RSQLite_2.3.4                 base64_2.0.1                  filelock_1.0.3                fansi_1.0.6                   httr_1.4.7                   
# [43] abind_1.4-5                   compiler_4.3.1                beanplot_1.3.1                rngtools_1.5.2                bit64_4.0.5                   BiocParallel_1.36.0           DBI_1.2.0                    
# [50] HDF5Array_1.30.0              biomaRt_2.58.0                MASS_7.3-60                   openssl_2.1.1                 rappdirs_0.3.3                DelayedArray_0.28.0           rjson_0.2.21                 
# [57] tools_4.3.1                   interactiveDisplayBase_1.40.0 httpuv_1.6.13                 glue_1.6.2                    quadprog_1.5-8                restfulr_0.0.15               promises_1.2.1               
# [64] nlme_3.1-164                  rhdf5filters_1.14.1           grid_4.3.1                    reshape2_1.4.4                generics_0.1.3                gtable_0.3.4                  tzdb_0.4.0                   
# [71] preprocessCore_1.64.0         tidyr_1.3.0                   hms_1.1.3                     xml2_1.3.6                    utf8_1.2.4                    BiocVersion_3.18.1            pillar_1.9.0                 
# [78] stringr_1.5.1                 limma_3.58.1                  later_1.3.2                   genefilter_1.84.0             splines_4.3.1                 dplyr_1.1.4                   lattice_0.22-5               
# [85] survival_3.5-7                rtracklayer_1.62.0            bit_4.0.5                     GEOquery_2.70.0               annotate_1.80.0               tidyselect_1.2.0              wheatmap_0.2.0               
# [92] scrime_1.3.5                  statmod_1.5.0                 stringi_1.8.3                 yaml_2.3.8                    codetools_0.2-19              tibble_3.2.1                  BiocManager_1.30.22          
# [99] cli_3.6.2                     xtable_1.8-4                  munsell_0.5.0                 Rcpp_1.0.11                   png_0.1-8                     XML_3.99-0.16                 ellipsis_0.3.2               
# [106] ggplot2_3.4.4                 readr_2.1.4                   blob_1.2.4                    prettyunits_1.2.0             mclust_6.0.1                  doRNG_1.8.6                   sparseMatrixStats_1.14.0     
# [113] bitops_1.0-7                  scales_1.3.0                  illuminaio_0.44.0             purrr_1.0.2                   crayon_1.5.2                  rlang_1.1.2                   KEGGREST_1.42.0           


