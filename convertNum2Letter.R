# this function is for converting the rack (8rowx12col) number nomenclature into letter annotation
# example: spot 13 = B1, spot 3 = A3, spot 50 = E2

test <- 12

convertNum2Letter <- function(test) {
  
  test.div12 <- ceiling(test/12) # since there are 12 columns, divide by 12 to determine the row number (round up)
  test.div12.char <- as.character(test.div12) # convert to character for lookup in dictionary
  
  # create a dictionary creating the references for letters in boxes
  # install.packages("hash")
  # library(hash)
  
  ## hash-2.2.6 provided by Decision Patterns
  plateref <- hash() 
  # set values
  plateref [["1"]] <- "A"
  plateref [["2"]] <- "B"
  plateref [["3"]] <- "C"
  plateref [["4"]] <- "D"
  plateref [["5"]] <- "E"
  plateref [["6"]] <- "F"
  plateref [["7"]] <- "G"
  plateref [["8"]] <- "H"
  
  letter.col <- plateref[[test.div12.char]] # grab the correct letter correcponding to the row
  letter.col 
  
  num.row <- test - (12*(test.div12-1)) # determine the column number by subtracting the original number from 12 x the number of full rows
  
  letterID <- paste(letter.col, num.row, sep = "")
  
  return(letterID)
  
}

testing <- convertNum2Letter(test)
testing
