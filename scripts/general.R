#CRAN
library(tidyverse)
library(glue)
#Bioconductor
library(Biostrings)
library(ShortRead)

#Reads in a Fastq file from the Dir
file2 <- "/Users/chriskent/phage_seq/seq_2/R2/"
input <- data.frame(file_dir = paste0(file2, "phagestockA_S10_L001_R2_001.fastq"),
                    output_name = c("unpan"))

counter <- function(x){
  raw_data2 <- readFastq(x$file_dir) %>% 
    sread() %>% 
    DNAStringSet() %>% 
    reverseComplement()

  # Removes which do not contain sequence just before VR
contains_vr <- vcountPattern(DNAString("TATTCTCACT"), raw_data2, max.mismatch = 1)
raw_data3 <- raw_data2[contains_vr == 1]

# Identifys the most common location of the sequence
vr_location <- vmatchPattern("TATTCTCACT", raw_data3) %>% unlist() 
vr_distro <- as.matrix(vr_location)[,1] %>% table()
max_distro <- names(vr_distro)[vr_distro == max(vr_distro)] %>%
  as.numeric()

# Generates the String to remove from left side, using the most common start + 11
l_str <- consensusString(Biostrings::subseq(raw_data3, start = 1, width = max_distro + 11))
test <- trimLRPatterns(Lpattern = l_str,
               subject = raw_data2,
               max.Lmismatch = 10)

# Identifys which sequences have been cut
seq_length <- table(width(test))
com_length <- names(seq_length)[seq_length == max(seq_length)] %>%
  as.numeric()

cut <- test[width(test) == com_length] %>% Biostrings::subseq(start = 1,
                                                      width = 36)

vr <- cut %>% translate()

tabler <- function(x){
  i <- table(x) %>% 
    as.data.frame() %>%
    arrange(by=desc(Freq))
  i <- DataFrame(seq = AAStringSet(i[,1]), 
                 freq = i[,2])
  i
}

count_dat <- tabler(vr)
write_csv(data.frame(count_dat),
          paste0("/Users/chriskent/phage_seq/count_dat/", x$output_name, "_autogen.csv"))

# Generates a log for the data
log <- glue::glue("For {x$output_name}
          Inital Length: {length(raw_data2)}
           Contains one VR: {length(raw_data3)}
           vr location: {max_distro}
           left consensus: {l_str}
           length: {length(test[width(test) == com_length])}
           {(length(test[width(test) == com_length])/length(raw_data2))*100}% sequences used")
write_file(log, 
           paste0("/Users/chriskent/phage_seq/count_dat/", x$output_name, "_log.csv"))

}
counter(input)
