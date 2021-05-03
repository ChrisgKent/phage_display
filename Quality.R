
barcodes_q <- raw_fastq@quality@quality %>% Biostrings::subseq(start = 5, width = 6) 
test <- PhredQuality(barcodes_q)
test <- as(test, "NumericList")
