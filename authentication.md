
library(Rsamtools)

library(magrittr)

library(ggplot2)

library(ggpubr)



BAM="C:/Users/Kursat/proje/İznik/Patojen_otantikasyonu/itk095-b1e1l1p1_TAATGCG-TTACTTA_L002_ARmerged.210803_A00187_0543_BHCFVTDRXY.all.fastq.gz_PipelineOutput/1302.output.sorted.bam"

mapdamage="C:/Users/Kursat/proje/İznik/Patojen_otantikasyonu/itk095-b1e1l1p1_TAATGCG-TTACTTA_L002_ARmerged.210803_A00187_0543_BHCFVTDRXY.all.fastq.gz_PipelineOutput/1302.output.bam/"

param <- ScanBamParam(tag = "NM", what = c("qwidth"))

bam <- scanBam(BAM, param = param)

edit <- bam[[1]]$tag$NM %>% table() %>% as.data.frame()
colnames(edit) <- c("NM", "Frequency")
edit$NM <- as.integer(as.character(edit$NM ))


edit <- edit[edit$NM <= 10,]


edit_plot <- ggplot(data = edit, aes(x=NM, y = Frequency)) + 
  geom_bar(stat="identity") + 
  xlab("Edit Distance")

length <- bam[[1]]$qwidth %>% table() %>% as.data.frame()

colnames(length) <- c("Length", "Frequency")
length$Length <- as.integer(as.character(length$Length))

length_plot <- ggplot(data = length, aes(x=Length, y = Frequency)) + 
  geom_bar(stat="identity") 


deamination5p <- read.table(paste0(mapdamage, "5pCtoT_freq.txt"), header=T, col.names = c("Position", "C_T"))
deamination3p <- read.table(paste0(mapdamage, "3pGtoA_freq.txt"), header=T, col.names = c("Position", "G_A"))

plot_5p <- ggplot(deamination5p,aes(x=Position, y=C_T)) + geom_line(color="red") +  xlab("Distance from 5' end") + ylab("C>T transition rate")
plot_3p <- ggplot(deamination3p,aes(x=Position, y=G_A)) + geom_line(color="blue") +  xlab("Distance from 3' end") + ylab("G>A transition rate") + scale_x_reverse() + 
  scale_y_continuous(position = "right")

last_plot <- ggarrange(edit_plot, length_plot, plot_5p, plot_3p, labels = c("a","b","c"))
