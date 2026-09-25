#!/usr/bin/env Rscript

library(ggplot2)
library(ggfortify)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
file.input <- args[1] # input _df
file.output <- args[2] # output name for fig in pdf format
my_sample <- args[3] # sample name (prefix) to color
# file.output2 <- args[3]
# file_rdata <- args[4]
df1=read.table(file.input,header=TRUE,sep=";",row.name=1, check.names = FALSE)



##delete col where the sum of the col is egal to 0, otherwise the PCA does not work
# data.without0val<-bin.bind[colSums(bin.bind) > 0]
data.without0val<-df1[colSums(df1) > 0]


##delete col where variance is equal to 0 otherwise it's not working for the PCA
data.without0val <- data.without0val[,apply(data.without0val, 2, var, na.rm=TRUE) != 0]
data.without0val$col_color <- ifelse(rownames(df1) == my_sample, "target", "untarget")
data.without0val <- data.without0val[order(data.without0val$col_color, decreasing = TRUE),]
data.without0val$col_color <- factor(data.without0val$col_color, levels = c("target","untarget"))



###
pca.1<-prcomp(data.without0val[,1:ncol(data.without0val)-1], center = TRUE, scale. = TRUE)
# pca.1<-prcomp(data.without0val, center = TRUE, scale. = TRUE)

# saveRDS(pca.1, file = file_rdata)



# pdf(file = file.output, width = 9, height = 9, pointsize = 10)
# autoplot(pca.1, data = data.without0val, colour = "col_color") +
#     scale_color_manual(values = c("red","black")) +
#     theme(legend.position = "none")
# dev.off
svg(file = file.output, width = 9, height = 9)
autoplot(pca.1, data = data.without0val, colour = "col_color", size = "col_color") +
    scale_color_manual(values = c("red","black")) +
    scale_size_manual(values = c(7,2)) +
    theme(legend.position = "none")
dev.off

# pdf(file = file.output2, width = 9, height = 9, pointsize = 10)
# autoplot(pca.1, data = data.without0val,label = TRUE, shape = FALSE)
# dev.off
###


