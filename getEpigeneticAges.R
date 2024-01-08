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
