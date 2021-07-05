#CRAN
library(tidyverse)
library(glue)
#Bioconductor
library(Biostrings)
library(ShortRead)

#Reads in a Fastq file from the Dir
file2 <- "/Users/chriskent/phage_seq/seq_2/R2/"
x <- data.frame(file_dir = paste0(file2, "phagestockA_S10_L001_R2_001.fastq"),
                    output_name = c("unpan"))

counter <- function(x){
  raw_data2 <- readFastq(x$file_dir) %>% 
    sread() %>% 
    DNAStringSet() %>% 
    reverseComplement()

# Removes which do not contain sequence just before VR
contains_vr <- vcountPattern(DNAString("TATTCTCACTCT"), raw_data2, max.mismatch = 1)
raw_data3 <- raw_data2[contains_vr == 1]

# Identifys the most common location of the sequence
vr_location <- vmatchPattern("TATTCTCACTCT", raw_data3) %>% unlist() 
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


# Recording results:
dir <- paste0("/Users/chriskent/phage_seq/count_dat/", x$output_name)
dir.create(dir)

write_csv(data.frame(count_dat),
          paste0(dir, "/", x$output_name, "_autogen.csv"))

write_csv(data.frame(seq = cut, aa = translate(cut)),
          paste0(dir, "/", "seq_", x$output_name, "_autogen.csv"))

# Generates a log for the data
log <- glue::glue("For {x$output_name}
          Inital Length: {length(raw_data2)}
           Contains one VR: {length(raw_data3)}
           vr location: {max_distro}
           left consensus: {l_str}
           length: {length(test[width(test) == com_length])}
           {(length(test[width(test) == com_length])/length(raw_data2))*100}% sequences used")
write_file(log, 
           paste0(dir, "/", x$output_name, "_log.csv"))



# plotting
test3 <- consensusMatrix(raw_data2)  %>% .[1:4,] %>%
  t() %>% as.data.frame() %>%
  mutate(sum = rowSums(.)) %>%
  mutate(position = 1:151) %>%
  select(position, 1:5)

test_p <- test3 %>% 
  pivot_longer(2:5, names_to = "bases", values_to = "count") %>%
  mutate(prop = count / sum)

plot_freq <- ggplot(test_p, aes(position,prop, col = bases))+
  geom_line()+
  theme_bw()+
  labs(title = paste0("Base Frequency of ", x$output_name),
       x = "Base position",
       y = "Proportion of each base")+
  scale_x_continuous(breaks = seq(0,200,10))+
  scale_color_discrete(name = "Bases")+
  theme(plot.title = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(colour = "grey"))

ggplot2::ggsave(filename = paste0(dir, "/", x$output_name, "_freq"),
                plot = plot_freq,
                device = "png",
                width = 10,
                height = 5,
                units = "in")
}


counter(input)


