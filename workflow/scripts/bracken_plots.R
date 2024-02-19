log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("ggplot2")
library("tidyr")

#cutoff
cutoff <- snakemake@params["abundance_cutoff"]

#Import data 
df <- read.table(snakemake@input[[1]], sep = "\t", header = T)

#Pivot longer

df.long <- df %>% pivot_longer(!name, names_to = "sample", values_to = "abund")

#For species less than the cutoff, change the name to other
df.long$name[df.long$abund<cutoff] <- "Other"


#Make the plot
plot<-
ggplot(subset(df.long, name!="Other"), aes(x=sample, y=abund, fill=name))+
  geom_bar(position = "stack", stat="identity")

ggsave(snakemake@output[[1]],plot)
ggsave(snakemake@output[[2]],plot)