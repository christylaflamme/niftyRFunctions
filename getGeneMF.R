# getGeneMF
# 02-24-22
# Christy LaFlamme
# Purpose: takes a list of genes and adds gene functional information 

library(mygene)
library(data.table)
library(readxl)

getGeneMF <- function(genelist, workDir, queryName) {
  
  # query functional information
  res <- queryMany(genelist, scopes = 'symbol', fields = c('entrezgene', 'go'), species = 'human')
  
  # create dataframe
  res.df <- data.frame(res) 
  res.pos <- res.df[is.na(res$notfound),] # only keep rows for which annotation information was found
  
  # write a loop that accumulates MF info for all successful queries
  for (i in 1:nrow(res.pos)){ # loop over the queries
    gene <- NULL
    res.MF <- NULL
    
    gene <- res.pos$query[i] # loop over the genes
    res.MF <- res[res$query == gene, 'go.MF'][[1]] # geet MF information
    res.pos$MF_summary[i] <- toString(res.MF$term) # add MF information to the dataframe
    
  }
  
  res.print <- res.pos[,c(1,9)] # print the gene name and MF information
  
  write.table(res.print, file = paste(workDir, queryName, sep = ""), quote = F, sep = "\t") # write to text file 
}
