############################################################################################
############################################################################################
# how to turn a dataframe into ExpansionHunter input catalog (R dataframe to json)
# LocusID
# LocusStructure
# ReferenceRegion
# VariantType

##############################
# install.packages("jsonlite")
library(jsonlite)

##############################
# imort dataframe
catalog.in <- as.data.frame(fread("./catalog.input.txt"))

# convert to json format
catalog.out <- toJSON(catalog.in)
write(catalog.out, file = "./catalog.output.json")

# repeat motifs need to be modified to include information about locus structure ()*

