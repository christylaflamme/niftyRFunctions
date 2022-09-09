# Running dmpFinder to find DMPs between condition and control towards deriving a methyation signature
# Christy LaFlamme
# Function created on 08-17-22

# load packages
library(data.table)
library(dplyr)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)
registerDoParallel(cores = 3)
library(pheatmap)

############################################################################################################
# Function to get number summaries of the derive and test cohorts

getCohortSummaries <- function(phenoData.derive, phenoData.test, condition.col, condition, outputDir, filename) { 
  
  # create a dataframe to store the summaries
  df <- data.frame(matrix(ncol = 2, nrow = 6))
  colnames(df) <- c("var", "num")
  var = c("derive_n", "test_n", "derive_n_condition", "derive_n_controls", "test_n_condition", "test_n_controls") # number of samples in the derive and test cohorts of the condition and controls
  df$var <- var
  
  # sum total
  derive_n = nrow(phenoData.derive) # total number of samples in the derive cohort
  test_n = nrow(phenoData.test) # total numbeer of samples in the test cohort
  print(paste0("There are ", as.character(derive_n), " total samples in the derive cohort and ", as.character(test_n), " total samples in the test cohort.")) 
  
  # derive
  derive_n_condition = length(which(phenoData.derive[[condition.col]] == condition)) # number of samples in the derive cohort that fall under the condition
  derive_n_controls = length(which(!phenoData.derive[[condition.col]] == condition)) + length(which(is.na(phenoData.derive[[condition.col]]))) # number of controls = not condition or NA
  print(paste0("There are ", as.character(derive_n_condition), " ", condition, " samples in the derive cohort and ", as.character(derive_n_controls), " controls."))
  
  # test 
  test_n_condition = length(which(phenoData.test[[condition.col]] == condition)) # number of samples in the test cohort that fall under the condition
  test_n_controls = length(which(!phenoData.test[[condition.col]] == condition)) + length(which(is.na(phenoData.test[[condition.col]]))) # number of controls = not condition or NA
  print(paste0("There are ", as.character(test_n_condition), " ", condition, " samples in the test cohort and ", as.character(test_n_controls), " controls."))
  
  # add the numerical information to the dataframe
  num <- c(derive_n, test_n, derive_n_condition, derive_n_controls, test_n_condition, test_n_controls)
  df$num <- num
  write.table(df, file = paste0(outputDir, "/", filename, "_cohort_n_summaries.txt"), quote = F, sep = "\t")
  
  return(df) # return the dataframe
  
}

############################################################################################################
# Function to get number summaries of the derive and test cohorts

getTestCohortSummaries <- function(phenoData.test, condition.col, condition, outputDir, filename) { 
  
  # create a dataframe to store the summaries
  df <- data.frame(matrix(ncol = 2, nrow = 3))
  colnames(df) <- c("var", "num")
  var = c("test_n", "test_n_condition", "test_n_controls") # number of samples in the test cohorts of the condition and controls
  df$var <- var
  
  # sum total
  derive_n = nrow(phenoData.derive) # total number of samples in the derive cohort
  test_n = nrow(phenoData.test) # total numbeer of samples in the test cohort
  print(paste0("There are ", as.character(test_n), " total samples in the test cohort.")) 
  
  # test 
  test_n_condition = length(which(phenoData.test[[condition.col]] == condition)) # number of samples in the test cohort that fall under the condition
  test_n_controls = length(which(!phenoData.test[[condition.col]] == condition)) + length(which(is.na(phenoData.test[[condition.col]]))) # number of controls = not condition or NA
  print(paste0("There are ", as.character(test_n_condition), " ", condition, " samples in the test cohort and ", as.character(test_n_controls), " controls."))
  
  # add the numerical information to the dataframe
  num <- c(test_n, test_n_condition, test_n_controls)
  df$num <- num
  write.table(df, file = paste0(outputDir, "/", filename, "_cohort_n_summaries.txt"), quote = F, sep = "\t")
  
  return(df) # return the dataframe
  
}

############################################################################################################
# Function to take an RGSet object, but the object for the samples of interest and get the beta values for analysis

getGRset <- function(RGSet, phenoData, manifest, filename, suffix) {
  
  print("Cutting the RGSset object.")
  # cut the RGSet object
  index <- as.numeric(which(colnames(RGSet) %in% phenoData$Sample_ID)) # get indices
  RGSet.index <- RGSet[,index]
  print(paste0("The total number of samples in the cut RGSet object are ", ncol(RGSet.index), "."))
  
  print("Converting the RGSset object to GRset.")
  # convert to GRset object
  GRset.funnorm <- preprocessFunnorm(RGSet.index)
  GRset.funnorm <- dropLociWithSnps(GRset.funnorm)
  print(paste0("The total number of samples in the GRset object are ", ncol(GRset.funnorm), "."))
  
  print("Removing sex chromosomes and mitochondrial DNA from GRset object...")
  # remove sex chr from analysis
  autosomes <- !(featureNames(GRset.funnorm) %in% manifest$probeID[manifest$CpG_chrm %in% c("chrX", "chrY", "chrM")]) 
  GRset.funnorm <- GRset.funnorm[autosomes,]
  print("Removed.")
  
  saveRDS(GRset.funnorm, file = paste(outputDir, filename, "_", ncol(RGSet.index), "n_", suffix, "_", ".GRset.funnorm.RDS", sep = ""))
  # GRset.funnorm <- readRDS(paste(outputDir, "GRset.funnorm.RDS", sep = ""))
  
  return(GRset.funnorm)
  
}

############################################################################################################
# Function to perform the dmpFinder DMP analysis
FindDMP <- function(GRset.funnorm, condition.col, condition, cohort.sum, outputDir, filename, manifest) { 
  
  # perform DMP analysis using dmpFinder
  print("Setting the conditions...")
  # set the conditions
  phenoData.GRset.funnorm <- pData(GRset.funnorm)
  Condition <- phenoData.GRset.funnorm[[condition.col]]    # create object for solved gene column
  Condition[!Condition == condition] <- "Control"           # designate controls
  Condition[is.na(Condition)] <- "Control"                  # designate controls
  print(paste0("There are ", length(Condition[Condition == condition]), " ", condition, " samples, and ", length(Condition[Condition == "Control"]), " controls."))
  
  print("Getting the beta values...")
  # get beta values
  beta <- getBeta(GRset.funnorm)
  # write.table(beta, file = paste(outputDir, "GRset.funnorm.beta.txt", sep = ""), quote = F, sep = "\t")
  print(paste0("The total number of samples in the beta value object are ", ncol(beta), "."))
  
  print("Performing dmpFinder DMP analysis.")
  dmp <- dmpFinder(beta, pheno = Condition, type = "categorical")
  
  print("Adding chromosomal information to DMP list...")
  dmp$probeID <- rownames(dmp)
  dmp <- merge(dmp, manifest[,c("probeID", "CpG_chrm", "CpG_beg", "CpG_end", "gene")], by.x = "probeID", by.y = "probeID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  
  print("Saving txt file with DMP info...")
  # because this "dmp" object only contains one matrix; the R object can be exported as a txt file that contains all the important information
  # saveRDS(dmp, file = paste0(outputDir, "/", filename, "_", cohort.sum[cohort.sum$var == "derive_n_condition",]$num, "n_", "control_", cohort.sum[cohort.sum$var == "derive_n_controls",]$num, "n_preprocessFunnorm_dmpFinder.RDS"))
  write.table(dmp, file = paste0(outputDir, "/", filename, "_", cohort.sum[cohort.sum$var == "derive_n_condition",]$num, "n_", "control_", cohort.sum[cohort.sum$var == "derive_n_controls",]$num, "n_preprocessFunnorm_dmpFinder.txt"), quote = F, sep = "\t", row.names = F)
  # dmp <- as.data.frame(read.table(paste0(outputDir, "/", filename, "_preprocessFunnorm_dmpFinder.txt", sep = "")))
  
  return(dmp)
  
}
############################################################################################################
# Function to plot heatmap of dmps for select samples

plotHeatmapDMPs <- function(GRset.funnorm, dmp, phenoData, cohort.sum, outputDir, filename, suffix) {
  
  #############################################
  print("Setting the conditions...")
  # get phenoData for GRset oboject
  phenoData.GRset.funnorm <- pData(GRset.funnorm)
  
  # set the conditions
  Condition <- phenoData.GRset.funnorm[[condition.col]]     # create object for solved gene column
  Condition[!Condition == condition] <- "Control"           # designate controls
  Condition[is.na(Condition)] <- "Control"                  # designate controls
  
  # set the sample IDs
  Sample_ID <- phenoData.GRset.funnorm[["Sample_ID"]] 
  
  # designate the groups for plotting (i.e. a key denoting which sample IDs correspond to which group
  group <- as.data.frame(Condition)
  rownames(group) <- Sample_ID # Sample_ID is the Sentrix ID + row and column numbers separated by an underscore
  print(group)
  #############################################
  
  #############################################
  print("Getting the beta values...")
  # get beta values
  beta <- getBeta(GRset.funnorm)
  beta.top1000dmps <- beta[which(rownames(beta) %in% dmp$probeID[1:1000]),] # create beta object with only top 1000 dmps
  #############################################
  
  #############################################
  print("Plotting the heatmap...")
  # plot heatmap 
  filenamePDF = paste0(filename, "_", cohort.sum[cohort.sum$var == paste0(suffix, "_n_condition"),]$num, "n_", "control_", cohort.sum[cohort.sum$var == paste0(suffix, "_n_controls"),]$num, "n_preprocessFunnorm_dmpFinder_top1000_methsig_", suffix)
  pdf(paste(outputDir, "/", filenamePDF, "_heatmap.pdf" , sep = ""))
  pheatmap.top1000dmps <- pheatmap(beta.top1000dmps, # matrix to be graphed as heatmap
                                   show_rownames=F, # remove rownames
                                   show_colnames = F, # remove colnames
                                   annotation_col = group) # add column annotation for condition vs. test
  
  dev.off()
  
  print("Getting sample order...")
  # get sample order
  col.order <- pheatmap.top1000dmps$tree_col$order 
  samples <- phenoData[match(colnames(beta), phenoData$Sample_ID),]
  pheatmap.ordered <- samples[col.order,] 
  filenameTXT = paste0(filename, "_", cohort.sum[cohort.sum$var == paste0(suffix, "_n_condition"),]$num, "n_", "control_", cohort.sum[cohort.sum$var == paste0(suffix, "_n_controls"),]$num, "n_preprocessFunnorm_dmpFinder_top1000_methsig_ranked_samples_", suffix)
  write.table(pheatmap.ordered, file = paste(outputDir, "/", filenameTXT, ".txt", sep = ""), quote = F, sep = "\t")
  
  print("Getting probe order")
  # get probe order
  row.order <- pheatmap.top1000dmps$tree_row$order
  probes <- beta.top1000dmps[row.order, col.order]
  filenameTXT.probes = paste0(filename, "_", cohort.sum[cohort.sum$var == paste0(suffix, "_n_condition"),]$num, "n_", "control_", cohort.sum[cohort.sum$var == paste0(suffix, "_n_controls"),]$num, "n_preprocessFunnorm_dmpFinder_top1000_methsig_ranked_probes_", suffix)
  write.table(probes, file = paste(outputDir, "/", filenameTXT.probes, ".txt", sep = ""), quote = F, sep = "\t")
  
}

############################################################################################################
# Function to derive and test methylation signature on pre-set cohorts

deriveTestMethSig <- function(phenoData.derive, phenoData.test, condition.col, condition, RGSet, annotDir, outputDir, filename) {
  
  #########################GET SIZES OF COHORTS##########################
  print("Loading the manifest file")
  # load the manifest file
  manifest <- fread(paste0(annotDir, "EPIC.hg19.manifest.tsv.gz"))
  
  #########################GET SIZES OF COHORTS##########################
  print("Creating final output directory...")
  outputDir <- paste0(outputDir, "/", filename, "_derive_", nrow(phenoData.derive), "n_test_", nrow(phenoData.test), "n_preprocessFunnorm_dmpFinder")
  print(outputDir)
  
  if (!file.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  print("Final output directory created...")
  
  print("Writing phenoData to output files...")
  write.table(phenoData.derive, file = paste(outputDir, "/phenoData.derive.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  write.table(phenoData.test, file = paste(outputDir, "/phenoData.test.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  #######################################################################
  
  #########################GET SIZES OF COHORTS##########################
  print("Assessing sizes of the cohorts...")
  cohort.sum <- getCohortSummaries(phenoData.derive, phenoData.test, condition.col, condition, outputDir, filename)
  print("Cohort size assessed...")
  #######################################################################
  
  #########################GET SIZES OF COHORTS##########################
  print("Getting GRset object for derive cohort...")
  GRset.derive <- getGRset(RGSet, phenoData.derive, manifest, filename, "derive")
  print("Completed...")
  #######################################################################
  
  ##################DMP ANLALYSIS ON DERIVE COHORT######################
  print("Performing DMP analysis on the derive cohort...")
  dmp <- FindDMP(GRset.derive, condition.col, condition, cohort.sum, outputDir, filename, manifest)
  print("DMP analysis complete...")
  #######################################################################
  
  ###################PLOT HEATMAP ON DERIVE COHORT#######################
  print("Plotting top 1000 DMPs on the derive cohort...")
  plotHeatmapDMPs(GRset.derive, dmp, phenoData.derive, cohort.sum, outputDir, filename, "derive")
  print("Plotted...")
  #######################################################################
  
  #########################GET SIZES OF COHORTS##########################
  print("Getting GRset object for test cohort...")
  GRset.test <- getGRset(RGSet, phenoData.test, manifest, filename, "test")
  print("Completed...")
  #######################################################################
  
  ###################PLOT HEATMAP ON TEST COHORT#########################
  print("Plotting top 1000 DMPs on the test cohort...")
  plotHeatmapDMPs(GRset.test, dmp, phenoData.test, cohort.sum, outputDir, filename, "test")
  print("Plotted...")
  #######################################################################
  
  print("Done. Successfully completed.")
  
}
############################################################################################################
# Function to derive and test methylation signature on pre-set cohorts

testMethSig <- function(phenoData.test, dmp, condition.col, condition, RGSet, annotDir, outputDir, filename) {
  
  #########################GET SIZES OF COHORTS##########################
  print("Loading the manifest file")
  # load the manifest file
  manifest <- fread(paste0(annotDir, "EPIC.hg19.manifest.tsv.gz"))
  
  #########################GET SIZES OF COHORTS##########################
  print("Creating final output directory...")
  outputDir <- paste0(outputDir, "/", filename, "_derive_", nrow(phenoData.derive), "n_test_", nrow(phenoData.test), "n_preprocessFunnorm_dmpFinder")
  print(outputDir)
  
  if (!file.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  print("Final output directory created...")
  
  print("Writing phenoData to output files...")
  write.table(phenoData.test, file = paste(outputDir, "/phenoData.test.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  #######################################################################
  
  #########################GET SIZES OF COHORTS##########################
  print("Assessing sizes of the cohorts...")
  cohort.sum <- getTestCohortSummaries(phenoData.test, condition.col, condition, outputDir, filename)
  print("Cohort size assessed...")
  #######################################################################
  
  #########################GET SIZES OF COHORTS##########################
  print("Getting GRset object for test cohort...")
  GRset.test <- getGRset(RGSet, phenoData.test, manifest, filename, "test")
  print("Completed...")
  #######################################################################
  
  ###################PLOT HEATMAP ON TEST COHORT#########################
  print("Plotting top 1000 DMPs on the test cohort...")
  plotHeatmapDMPs(GRset.test, dmp, phenoData.test, cohort.sum, outputDir, filename, "test")
  print("Plotted...")
  #######################################################################
  
  print("Done. Successfully completed.")
  
}

############################################################################################################
# Example usage of the function(s)
############################################################################################################
# set number of samples
n = "1253"

# set working directories
methDir = "/home/claflamm/methylation/"
workDir = paste0("/home/claflamm/methylation/raw_idats/", n, "_all_metharray/output/", n, "_all_metharray_dev/", n, "_all_metharray_dev_preprocessIllumina/")
outputDir = "/home/claflamm/methylation2.0/SCN1A/dmpFinder/subset_SCN1A_vs_controls/"
SCN1ADir = "/home/claflamm/methylation2.0/SCN1A/"
annotDir <- "/home/claflamm/methylation/annotations/" # This directory should contain the file: EPIC.hg19.manifest.tsv.gz
# wget https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz

############################################################################################################
# load RGSet object and phenotype data
phenoData <- fread(paste0(methDir, "meta_data/mefford_methylation_pl1-11_and_controls.csv"), skip =7)
RGSet <- readRDS(paste(workDir, "raw_data/RGSet.all.RDS", sep = ""))

# filter for blood samples
phenoData.blood <- phenoData[phenoData$Sample_Type == "Blood" | phenoData$Sample_Type == "Unknown",]

# subset out the samples by gene
CHD2.samples <- phenoData.blood %>% filter(Solved_Gene == "CHD2") # all SCN1A samples - 49 samples
SCN1A.samples <- phenoData.blood %>% filter(Solved_Gene == "SCN1A") # all CHD2 samples - 16 samples

exclude_controls <- fread(paste0(outputDir, "1230_blood_saliva_PBMC_samples_0.995.0.005_sample.summary.txt"))
exclude_controls <- exclude_controls %>% filter(nDMRsperSample > 10) %>% select(Sample_ID)

# ############################################################################################################
# SCN1A NONSENSE
SCN1A.samples.nonsense <- SCN1A.samples %>% filter(Solved_Mutation_Type == "Nonsense")
# subset samples for nonsense group; starting with just plate 8 samples to mitigate batch effects
SCN1A.samples.derive <- SCN1A.samples.nonsense[3:12,] # 10 samples
SCN1A.samples.test <- SCN1A.samples.nonsense[13:16,] # 4 samples

# subset control samples that were run at St. Jude to help mitigate batch effects
control.samples <- phenoData.blood %>% filter(Sample_Plate == "mefford_sjcrh_dee_methylation_4~8") # | Sample_Plate == "mefford_sjcrh_dee_methylation_5~9" | Sample_Plate == "mefford_sjcrh_dee_methylation_6~10"
dim(control.samples)
control.samples <- control.samples %>% filter(!Solved_Gene == "SCN1A" | is.na(Solved_Gene)) # get rid of SCN1A samples
dim(control.samples)
control.samples <- control.samples %>% filter(Sample_Designation == "Proband_Unsolved" | Sample_Designation == "Proband_Solved" ) # get rid of parents etc.
dim(control.samples)
control.samples <- control.samples %>% filter(!Notes == "Unsolved Dravet" & !Notes == "DMRs_3000" | is.na(Notes)) # remove unsolved Dravet - these may have "missing" SCN1A variants; remove samples with a lot of DMRs
dim(control.samples)
control.samples <- control.samples %>% filter(!Sample_ID %in% exclude_controls)
dim(control.samples)

# nonsense derive 10 test 4
control.samples.derive <- control.samples[1:40,] # 40 samples
control.samples.test <- control.samples[41:52,]# 12 samples

phenoData.derive <- NA
phenoData.test <- NA

phenoData.derive <- rbind(SCN1A.samples.derive, control.samples.derive)
phenoData.test <- rbind(SCN1A.samples.test, control.samples.test)

# run function
condition.col = NA
condition = NA
filename = NA

condition.col = "Solved_Gene"
condition = "SCN1A"
filename = "SCN1A_nonsense"

deriveTestMethSig(phenoData.derive, phenoData.test, condition.col, condition, RGSet, annotDir, outputDir, filename)


############################################################################################################
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12                                     doParallel_1.0.17                                   IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
# [4] IlluminaHumanMethylationEPICmanifest_0.3.0          minfi_1.42.0                                        bumphunter_1.38.0                                  
# [7] locfit_1.5-9.6                                      iterators_1.0.14                                    foreach_1.5.2                                      
# [10] Biostrings_2.64.0                                   XVector_0.36.0                                      SummarizedExperiment_1.26.1                        
# [13] Biobase_2.56.0                                      MatrixGenerics_1.8.1                                matrixStats_0.62.0                                 
# [16] GenomicRanges_1.48.0                                GenomeInfoDb_1.32.2                                 IRanges_2.30.0                                     
# [19] S4Vectors_0.34.0                                    BiocGenerics_0.42.0                                 dplyr_1.0.9                                        
# [22] data.table_1.14.2                                  
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-3          rjson_0.2.21              ellipsis_0.3.2            siggenes_1.70.0           mclust_5.4.10             base64_2.0               
# [7] rstudioapi_0.13           bit64_4.0.5               AnnotationDbi_1.58.0      fansi_1.0.3               xml2_1.3.3                codetools_0.2-18         
# [13] splines_4.2.0             sparseMatrixStats_1.8.0   cachem_1.0.6              scrime_1.3.5              Rsamtools_2.12.0          annotate_1.74.0          
# [19] dbplyr_2.2.1              png_0.1-7                 HDF5Array_1.24.1          BiocManager_1.30.18       readr_2.1.2               compiler_4.2.0           
# [25] httr_1.4.3                assertthat_0.2.1          Matrix_1.4-1              fastmap_1.1.0             limma_3.52.2              cli_3.3.0                
# [31] prettyunits_1.1.1         tools_4.2.0               gtable_0.3.0              glue_1.6.2                GenomeInfoDbData_1.2.8    rappdirs_0.3.3           
# [37] doRNG_1.8.2               Rcpp_1.0.9                vctrs_0.4.1               rhdf5filters_1.8.0        multtest_2.52.0           preprocessCore_1.58.0    
# [43] nlme_3.1-157              rtracklayer_1.56.1        DelayedMatrixStats_1.18.0 stringr_1.4.0             lifecycle_1.0.1           restfulr_0.0.15          
# [49] rngtools_1.5.2            XML_3.99-0.10             beanplot_1.3.1            scales_1.2.0              zlibbioc_1.42.0           MASS_7.3-57              
# [55] hms_1.1.1                 rhdf5_2.40.0              GEOquery_2.64.2           RColorBrewer_1.1-3        yaml_2.3.5                curl_4.3.2               
# [61] memoise_2.0.1             biomaRt_2.52.0            reshape_0.8.9             stringi_1.7.8             RSQLite_2.2.14            genefilter_1.78.0        
# [67] BiocIO_1.6.0              GenomicFeatures_1.48.3    filelock_1.0.2            BiocParallel_1.30.3       rlang_1.0.3               pkgconfig_2.0.3          
# [73] bitops_1.0-7              nor1mix_1.3-0             lattice_0.20-45           purrr_0.3.4               Rhdf5lib_1.18.2           GenomicAlignments_1.32.0 
# [79] bit_4.0.4                 tidyselect_1.1.2          plyr_1.8.7                magrittr_2.0.3            R6_2.5.1                  generics_0.1.3           
# [85] DelayedArray_0.22.0       DBI_1.1.3                 pillar_1.7.0              survival_3.2-13           KEGGREST_1.36.2           RCurl_1.98-1.7           
# [91] tibble_3.1.7              crayon_1.5.1              utf8_1.2.2                BiocFileCache_2.4.0       tzdb_0.3.0                progress_1.2.2           
# [97] grid_4.2.0                blob_1.2.3                digest_0.6.29             xtable_1.8-4              tidyr_1.2.0               illuminaio_0.38.0        
# [103] munsell_0.5.0             openssl_2.0.2             askpass_1.1               quadprog_1.5-8           

