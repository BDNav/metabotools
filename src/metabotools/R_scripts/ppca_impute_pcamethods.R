

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
message(args[1])
col <- as.integer(args[2])   
message(args[2])

library(pcaMethods)
data <- read.csv(args[1])
data2 <- data[,-c(1:col)] # remove first X columns (group and name)
pc <- pca(data2, nPcs=3, method="ppca")
imputed <- completeObs(pc)
write.csv(imputed, "ppca_imputed.csv")