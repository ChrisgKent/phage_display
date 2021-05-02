#Run this file 

#This take a set of XStringSet, and returns the number times the sequnce appears.
counter <- function(xstring, mismatch = 0){
  unique_x <- unique(xstring)
  x <- vcountPDict(unique_x,xstring, collapse = 1, max.mismatch = mismatch)
  dat <- DataFrame(seq = unique_x, count = x)
  dat
}

subseter <- function(bc){
  bc <- DNAString(bc)
  c <- vcountPattern(bc, barcodes)
  raw_dat[c,]
}

