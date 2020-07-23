library(veriNA3d)
library(dplyr)
library(tidyverse)

source("full_data.RData")

#Extract all those which are TRUE for motifs
full_data2 <- full_data[which( full_data$motifs %in% TRUE),]

#TRUE - helix and non-helix - FALSE - all the rest

#Exctract from those the RNA associated ones (remove any D*)
full_data3 <- full_data2[grep("^[AUCG]-[AUGC]-[AUGC]$", full_data2$ntinfo.localenv, perl=T),]

#Create a factor table with all the trinucleotides and the ocurrences within the dataset
table <- table(as.factor(full_data3$ntinfo.localenv))
data <- as.data.frame(table)


#Create the data frame with ascending order data
sorted <- data.frame(reorder(data$Var1, data$Freq))
names(sorted) <- paste("Trinucleotide")
sorted_data <- data.frame (Nucleotide = sorted, Freq = data$Freq)

#Crete the plot with ascending order data
ggplot(sorted_data, aes(x = Trinucleotide, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90,margin = margin(1, unit = "cm"),vjust =1))
