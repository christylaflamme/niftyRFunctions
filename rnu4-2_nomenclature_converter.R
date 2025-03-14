####################################################################################
# functions to convert g. and n. nomenclature of RNU4-2 variants (in hg38)
# Christy LaFlamme
# 03/14/25
# g=genomic
# n=nucleotide (RNA)
##############################################
# load libraries
# install.packages("hash")
library(hash)
##############################################
# RNU4-2
# from g --> n (SNVs)

coord <- "12:120291830:G:A"

g2n.snv.rnu4.2 <- function(coord="12:120291826:T:C") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if chromosome listed is chr12
  if (coord.split[[1]][1] == "12") {
    
    # calculate the position in the RNA (n.)
    n.coord <- 120291903 - as.numeric(coord.split[[1]][2]) + 1
    
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
    print("RNU4-2 is on chromosome 12. Please enter valid coordinates.")
    
  }
  
}


test <- g2n.snv.rnu4.2(coord)
##############################################
# RNU4-2
# from n --> g (SNVs)

coord <- "n:64:G:C"

n2g.snv.rnu4.2 <- function(coord="n:78:A:G") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if first split character is n
  if (coord.split[[1]][1] == "n") {
    
    # calculate the position in the RNA (n.)
    g.coord <- 120291903 - as.numeric(coord.split[[1]][2]) + 1
    
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
    return(paste("12:", g.coord, ":", g.ref, ":", g.alt, sep=""))
    
} else {
  print("Expecting n. nomenclature. Please enter valid coordinates.")

  }

}

test <- n2g.snv.rnu4.2(coord)

##############################################
# RNU4-2
# from g --> n (Indels)

coord <- "12:120291826:T:T:A"

g2n.indel.rnu4.2 <- function(coord="12:120291839:T:T:A") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if chromosome listed is chr12
  if (coord.split[[1]][1] == "12") {
    
    # calculate the position in the RNA (n.)
    n.coord.start <- 120291903 - as.numeric(coord.split[[1]][2])
    n.coord.end <- 120291903 - as.numeric(coord.split[[1]][2]) + 1
    
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
    print("RNU4-2 is on chromosome 12. Please enter valid coordinates.")

  }

}


test <- g2n.indel.rnu4.2(coord)

##############################################
# RNU4-2
# from n --> g (Indels)

coord <- "n:73:74:ins:A"

n2g.indel.rnu4.2 <- function(coord="n:64:65:ins:T") { # minimal error handling; assuming correct formatting is used!
  
  # split the coordinate by :
  coord.split <- strsplit(coord, ":")
  
  # check if first split character is n
  if (coord.split[[1]][1] == "n") {
    
    # calculate the position in the RNA (n.)
    g.coord <- 120291903 - as.numeric(coord.split[[1]][3]) + 1
    
    # set values for complementary bases (essentially a dictionary)
    baseref <- hash() 
    baseref [["A"]] <- "T"
    baseref [["T"]] <- "A"
    baseref [["C"]] <- "G"
    baseref [["G"]] <- "C"
    
    # index over the RNU4-2 sequence to find the REF (since n. nomenclature for indels doesn't include the REF)
    seq = "AGCTTTGCGCAGTGGCAGTATCGTAGCCAATGAGGTTTATCCGAGGCGCGATTATTGCTAATTGAAAACTTTTCCCAATACCCCGCCATGACGACTTGAAATATAGTCGGCATTGGCAATTTTTGACAGTCTCTACGGAGA"
    seq <- unlist(strsplit(seq, ""))
    
    # pull the complementary bases for REF and ALT given
    g.ref <- baseref[[seq[as.numeric(coord.split[[1]][3])]]]
    g.alt <- baseref[[coord.split[[1]][5]]]
    
    # return the converted n. format
    return(paste("12:", g.coord, ":", g.ref, ":", g.ref, g.alt, sep=""))
    
  } else {
    print("Expecting n. nomenclature. Please enter valid coordinates.")

  }

}

test <- n2g.indel.rnu4.2(coord)
