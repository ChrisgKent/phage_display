#CRAN
library(tidyverse)
#Bioconductor
library(Biostrings)
library(ShortRead)

#Reads in a Fastq file from the Dir
file <- "PRIVATE_DO_NOT_COMMIT/data/VIRSCAN1_S1_L001_R1_001.fastq"
raw_fastq <- readFastq(file)
#Extracts the sequence data from the file, as asigns it to a DNAStringSet
raw_dat <- sread(raw_fastq) %>% DNAStringSet()

#Extracts the section of the sequene where the barcode is expected.
#Counter determins the number of time a single barcode pattern has been read
barcodes <- Biostrings::subseq(raw_dat, start = 50, width = 6)
barcode_count <- counter(barcodes)


write_delim(as.data.frame(barcode_count), "PRIVATE_DO_NOT_COMMIT/files/barcode_count.txt", col_names = FALSE)

# This picks the 20 barcodes with the highest count, and assigns them to a list
# As strings for memory reasons.
top_20 <- barcode_count[1:200,]
seq_str <- top_20$seq %>% toString() %>% str_split(", ", simplify = TRUE)
blist <- as.list(seq_str)
names(blist) <- seq_str


# This applys the barcode list to the raw_data file, and subsets sequences based on the barcodes
bc_subset <- mclapply(blist, subseter) %>% DNAStringSetList()
