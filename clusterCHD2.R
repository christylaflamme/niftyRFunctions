###############################################################################

# clusterCHD2
# Christy LaFlamme 
# 02-25-22
# Purpose: cluster methylation array samples using CHD2 methylation signature probes against known CHD2 patients and controls 

###############################################################################
# INPUT
# phenoData.all: all of the phenotype data
  # Required columns:
  # Sample_Name
  # Sample_Type
  # Sample_ID (SentrixID)
# phenoData.test: phenoData information for the samples to be tested
# mVal: object containing the normalized M-values from the methylation array
# workDir: path to working directory
# filename: filename to be used for output files (no extensions)

# OUTPUT
# heatmap of clustering for CHD2 probes containing annotations for CHD2 samples, controls, and test samples
# text file of ranked samples with phenotype data

###############################################################################
# load dependencies
library(dplyr)
library(tidyverse)
library(minfi) 
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
library(pheatmap)
library(RColorBrewer)

annoEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annoEPIC)

# load CHD2 probelist
CHD2Dir = "/Users/claflamm/OneDrive - St. Jude Children's Research Hospital/Mefford/methylation/CHD2/"
CHD2.meth.sig.probes <- readRDS(file = paste(CHD2Dir, "CHD2.meth.sig.probes.RDS", sep = ""))

###############################################################################

clusterCHD2 <- function(phenoData.all, phenoData.test, mVal, workDir, filename) {
  
  # subset CHD2 patients
  phenoData.CHD2 <- phenoData %>% 
    filter(Sample_Type == "Blood" | Sample_Type == "Unknown") %>%
    filter(Solved_Gene == "CHD2")
  
  # subset controls (make sure not test of CHD2)
  phenoData.controls <- phenoData %>% 
    filter(!Sample_ID %in% phenoData.test$Sample_ID) %>%
    filter(!Sample_ID %in% phenoData.CHD2$Sample_ID) %>%
    filter(Sample_Type == "Blood" | Sample_Type == "Unknown") %>%
    filter(!Solved_Gene == "CHD2" | is.na(Solved_Gene))
  
  # select 150 blood samples randomly as controls
  phenoData.controls <- phenoData.controls[sample(1:nrow(phenoData.controls), 150), ] 
  
  # combine all sample data
  phenoData.combined <- rbind(phenoData.CHD2, phenoData.test, phenoData.controls)
  
  # get mVal for selected samples and CHD2 methylation signature probes
  mVal <- mVal[rownames(mVal) %in% CHD2.meth.sig.probes, which(colnames(mVal) %in% phenoData.combined$Sample_ID)]
  
  # remove Na/NaN/Inf
  mVal<- mVal[is.finite(rowSums(mVal)),]
  
  # set groups for annotation
  group <- as.data.frame(phenoData.combined$Solved_Gene)
  colnames(group)[1] <- "Group"
  row.names(group) <- phenoData.combined$Sample_ID
  group[rownames(group) %in% phenoData.CHD2$Sample_ID, "Group"] <- "CHD2" 
  group[rownames(group) %in% phenoData.test$Sample_ID, "Group"] <- "Test" 
  group[rownames(group) %in% phenoData.controls$Sample_ID, "Group"] <- "Control" 
  
  # plot heatmap
  pdf(paste(workDir, filename, "_heatmap.pdf" , sep = ""))
  pheatmap.CHD2 <- pheatmap(mVal, # matrix to be graphed as heatmap
                            show_rownames=F, # remove rownames 
                            show_colnames = F, # rwmove colnames
                            annotation_col = group) # add column annotation for CHD2 vs. test
  
  dev.off()
  
  # get list of samples in the order of the mVal columns
  samples <- phenoData.combined[match(colnames(mVal), phenoData.combined$Sample_ID),] 
  
  # get ranked list of samples based on the hierarchical clustering
  col.order <- pheatmap.CHD2$tree_col$order # column order of the heatmap
  pheatmap.ordered <- samples[col.order,] # order to rows of the samples
  write.table(pheatmap.ordered, file = paste(workDir, filename, "_ranked_samples.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  
}

###############################################################################
# example usage

phenoData.all <- as.data.frame(read.metharray.sheet(dataDir))
phenoData.test <- phenoData.all %>%
  filter(Sample_Plate == "2632-GC-P1")
mVal <- readRDS(paste(outputDir, "933_all_DEE_controls_blood_saliva/mVal.filtered.RDS", sep = "")) 
workDir = "/Users/claflamm/OneDrive - St. Jude Children's Research Hospital/Mefford/methylation/CHD2/output/Gemma_Fb_iPSC_NPC/"
filename = "clusterCHD2_Gemma_Fb_iPSC_NPC"
clusterCHD2 <- clusterCHD2(phenoData.all, phenoData.test, mVal, workDir, filename)

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
#   [1] pheatmap_1.0.12                                     Rtsne_0.15                                          limma_3.48.3                                       
# [4] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0 IlluminaHumanMethylationEPICmanifest_0.3.0          minfi_1.38.0                                       
# [7] bumphunter_1.34.0                                   locfit_1.5-9.4                                      iterators_1.0.13                                   
# [10] foreach_1.5.1                                       Biostrings_2.60.2                                   XVector_0.32.0                                     
# [13] SummarizedExperiment_1.22.0                         MatrixGenerics_1.4.3                                matrixStats_0.61.0                                 
# [16] forcats_0.5.1                                       stringr_1.4.0                                       purrr_0.3.4                                        
# [19] readr_2.1.0                                         tidyr_1.1.4                                         tibble_3.1.6                                       
# [22] ggplot2_3.3.5                                       tidyverse_1.3.1                                     dplyr_1.0.7                                        
# [25] readxl_1.3.1                                        data.table_1.14.2                                   mygene_1.28.0                                      
# [28] GenomicFeatures_1.44.2                              AnnotationDbi_1.54.1                                Biobase_2.52.0                                     
# [31] GenomicRanges_1.44.0                                GenomeInfoDb_1.28.4                                 IRanges_2.26.0                                     
# [34] S4Vectors_0.30.2                                    BiocGenerics_0.38.0                                
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.3.0           Hmisc_4.6-0               BiocFileCache_2.0.0       plyr_1.8.6                splines_4.1.1             BiocParallel_1.26.2      
# [7] digest_0.6.28             htmltools_0.5.2           fansi_0.5.0               magrittr_2.0.1            checkmate_2.0.0           memoise_2.0.0            
# [13] cluster_2.1.2             tzdb_0.2.0                annotate_1.70.0           modelr_0.1.8              askpass_1.1               siggenes_1.66.0          
# [19] prettyunits_1.1.1         jpeg_0.1-9                colorspace_2.0-2          blob_1.2.2                rvest_1.0.2               rappdirs_0.3.3           
# [25] haven_2.4.3               xfun_0.28                 crayon_1.4.2              RCurl_1.98-1.5            jsonlite_1.7.2            genefilter_1.74.1        
# [31] GEOquery_2.60.0           survival_3.2-13           glue_1.5.0                gtable_0.3.0              zlibbioc_1.38.0           DelayedArray_0.18.0      
# [37] Rhdf5lib_1.14.2           HDF5Array_1.20.0          scales_1.1.1              DBI_1.1.1                 rngtools_1.5.2            Rcpp_1.0.7               
# [43] xtable_1.8-4              progress_1.2.2            htmlTable_2.3.0           foreign_0.8-81            bit_4.0.4                 mclust_5.4.8             
# [49] preprocessCore_1.54.0     Formula_1.2-4             sqldf_0.4-11              htmlwidgets_1.5.4         httr_1.4.2                RColorBrewer_1.1-2       
# [55] ellipsis_0.3.2            farver_2.1.0              pkgconfig_2.0.3           reshape_0.8.8             XML_3.99-0.8              nnet_7.3-16              
# [61] dbplyr_2.1.1              utf8_1.2.2                tidyselect_1.1.1          rlang_0.4.12              munsell_0.5.0             cellranger_1.1.0         
# [67] tools_4.1.1               cachem_1.0.6              cli_3.1.0                 gsubfn_0.7                generics_0.1.1            RSQLite_2.2.8            
# [73] broom_0.7.10              fastmap_1.1.0             yaml_2.2.1                knitr_1.36                bit64_4.0.5               fs_1.5.0                 
# [79] beanplot_1.2              scrime_1.3.5              KEGGREST_1.32.0           nlme_3.1-153              doRNG_1.8.2               sparseMatrixStats_1.4.2  
# [85] nor1mix_1.3-0             xml2_1.3.2                biomaRt_2.48.3            compiler_4.1.1            rstudioapi_0.13           filelock_1.0.2           
# [91] curl_4.3.2                png_0.1-7                 reprex_2.0.1              stringi_1.7.5             lattice_0.20-45           Matrix_1.3-4             
# [97] multtest_2.48.0           vctrs_0.3.8               pillar_1.6.4              lifecycle_1.0.1           rhdf5filters_1.4.0        bitops_1.0-7             
# [103] rtracklayer_1.52.1        R6_2.5.1                  BiocIO_1.2.0              latticeExtra_0.6-29       gridExtra_2.3             codetools_0.2-18         
# [109] MASS_7.3-54               assertthat_0.2.1          chron_2.3-56              rhdf5_2.36.0              proto_1.0.0               openssl_1.4.5            
# [115] rjson_0.2.20              withr_2.4.2               GenomicAlignments_1.28.0  Rsamtools_2.8.0           GenomeInfoDbData_1.2.6    hms_1.1.1                
# [121] quadprog_1.5-8            grid_4.1.1                rpart_4.1-15              base64_2.0                DelayedMatrixStats_1.14.3 illuminaio_0.34.0        
# [127] lubridate_1.8.0           base64enc_0.1-3           restfulr_0.0.13     
