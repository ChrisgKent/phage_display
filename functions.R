#Run this file 

#This take a set of XStringSet, and returns the number times the sequnce appears.
# mc arg. allows the use of mutlicore processing. Does throw a warning if the unique(xstring) isnt devisable by 10
counter <- function(xstring, mismatch = 0, mc = FALSE){
  unique_x <- unique(xstring)
  if(mc == TRUE){
    list <- split(unique_x, 1:10)
    x <- mclapply(list, function(x){
      vcountPDict(x,xstring, collapse = 1, max.mismatch = mismatch)})
    dat <- DataFrame(seq = unlist(list), count = unlist(x))
    
  }else{x <- vcountPDict(unique_x,xstring, collapse = 1, max.mismatch = mismatch)
  dat <- DataFrame(seq = unique_x, count = x)
  }
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


