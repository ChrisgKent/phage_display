#Run this file 

#This take a set of XStringSet, and returns the number times the sequnce appears.
counter <- function(xstring, mismatch = 0){
  unique_x <- unique(xstring)
  x <- vcountPDict(unique_x,xstring, collapse = 1, max.mismatch = mismatch)
  dat <- DataFrame(seq = unique_x, count = x)
  rev(dat[(order(dat$count)),])
}


subseter <- function(bc){
  bc <- DNAString(bc)
  c <- countPDict(barcodes, bc, max.mismatch = 0) %>% as.logical()
  raw_dat[c,]
}

vr_extract <- function(x){
  i <- x %>% Biostrings::subseq(start= 57, width = 36)
  Biostrings::translate(i)
}

test <- mclapply(bc_subset, vr_extract)
test2 <- mclapply(test, counter)
test2
