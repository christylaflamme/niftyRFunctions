####################################################################################
# functions to convert g. and n. nomenclature of RNU2-2(P) variants (in hg38)
# Christy LaFlamme
# 03/14/25
# g=genomic
# n=nucleotide (RNA)
##############################################
# load libraries
# install.packages("hash")
library(hash)
##############################################
# RNU2-2
# from g --> n (SNVs)

coord <- "11:62841806:C:T"

g2n.snv.rnu2.2 <- function(coord="11:62841810:A:T") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if chromosome listed is chr11
  if (coord.split[[1]][1] == "11") {
    
    # calculate the position in the RNA (n.)
    n.coord <- 62841809 - as.numeric(coord.split[[1]][2]) + 1
    
    # set values for complementary bases (essentially a dictionary)
    baseref <- hash() 
    baseref [["A"]] <- "T"
    baseref [["T"]] <- "A"
    baseref [["C"]] <- "G"
    baseref [["G"]] <- "C"
    
    # pull the complementary bases for REF and ALT given
    n.ref <- baseref [[coord.split[[1]][3]]]
    n.alt <- baseref [[coord.split[[1]][4]]]
    
    # return the converted n. format
    return(paste("n.", n.coord, n.ref, ">", n.alt, sep=""))
    
  } else {
    print("RNU2-2 is on chromosome 11. Please enter valid coordinates.")

  }

}


test <- g2n.snv.rnu2.2(coord)
##############################################
# RNU2-2
# from n --> g (SNVs)

coord <- "n:4:G:A"

n2g.snv.rnu2.2 <- function(coord="n:4:G:A") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if first split character is n
  if (coord.split[[1]][1] == "n") {
    
    # calculate the position in the RNA (n.)
    g.coord <- 62841809 - as.numeric(coord.split[[1]][2]) + 1
    
    # set values for complementary bases (essentially a dictionary)
    baseref <- hash() 
    baseref [["A"]] <- "T"
    baseref [["T"]] <- "A"
    baseref [["C"]] <- "G"
    baseref [["G"]] <- "C"
    
    # pull the complementary bases for REF and ALT given
    g.ref <- baseref [[coord.split[[1]][3]]]
    g.alt <- baseref [[coord.split[[1]][4]]]
    
    # return the converted n. format
    return(paste("11:", g.coord, ":", g.ref, ":", g.alt, sep=""))
    
  } else {
    print("Expecting n. nomenclature. Please enter valid coordinates.")
    
  }
  
}

test <- n2g.snv.rnu2.2(coord)

##############################################
# RNU2-2
# from g --> n (Indels)

coord <- "11:62841806:C:C:T"

g2n.indel.rnu2.2 <- function(coord="11:62841806:C:C:T") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if chromosome listed is chr11
  if (coord.split[[1]][1] == "11") {
    
    # calculate the position in the RNA (n.)
    n.coord.start <- 62841809 - as.numeric(coord.split[[1]][2])
    n.coord.end <- 62841809 - as.numeric(coord.split[[1]][2]) + 1
    
    # set values for complementary bases (essentially a dictionary)
    baseref <- hash() 
    baseref [["A"]] <- "T"
    baseref [["T"]] <- "A"
    baseref [["C"]] <- "G"
    baseref [["G"]] <- "C"
    
    # pull the complementary bases for REF and ALT given
    n.alt <- baseref [[coord.split[[1]][5]]]
    
    # return the converted n. format
    return(paste("n.", n.coord.start, "_", n.coord.end, "ins", n.alt, sep=""))
    
  } else {
    print("RNU2-2P is on chromosome 11. Please enter valid coordinates.")
    
  }
  
}


test <- g2n.indel.rnu2.2(coord)

##############################################
# RNU2-2
# from n --> g (Indels)

coord <- "n:3:4:ins:A"

n2g.indel.rnu2.2 <- function(coord="n:3:4:ins:A") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if first split character is n
  if (coord.split[[1]][1] == "n") {
    
    # calculate the position in the RNA (n.)
    g.coord <- 62841809 - as.numeric(coord.split[[1]][3]) + 1
    
    # set values for complementary bases (essentially a dictionary)
    baseref <- hash() 
    baseref [["A"]] <- "T"
    baseref [["T"]] <- "A"
    baseref [["C"]] <- "G"
    baseref [["G"]] <- "C"
    
    # index over the RNU2-2 sequence to find the REF (since n. nomenclature for indels doesn't include the REF)
    seq = "ATCGCTTCTCGGCCTTTTGGCTAAGATCAAGTGTAGTATCTGTTCTTATCAGTTTAATATCTGATACGTCCTCTATCCGAGGACAATATATTAAATGGATTTTTGGAAATAGGAGATGGAATAGGAGCTTGCTCCGTCCACTCCACGCATCGACCTGGTATTGCAGTACTTCCAGGAACGGTGCACTCTCC"
    seq <- unlist(strsplit(seq, ""))
    
    # pull the complementary bases for REF and ALT given
    g.ref <- baseref[[seq[as.numeric(coord.split[[1]][3])]]]
    g.alt <- baseref[[coord.split[[1]][5]]]
    
    # return the converted n. format
    return(paste("11:", g.coord, ":", g.ref, ":", g.ref, g.alt, sep=""))
    
  } else {
    print("Expecting n. nomenclature. Please enter valid coordinates.")
    
  }
  
}

test <- n2g.indel.rnu2.2(coord)
