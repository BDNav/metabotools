

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
message(args)   

library("gWQS")

data <- read.table(file=args[1], header=TRUE, sep=',')


PCBs<-names(data)[3:ncol(data)]

results <- gwqs(Group ~ wqs, mix_name = PCBs, data = data, q = 10, validation = 0.2, b = strtoi(args[2]), b1_pos = FALSE, b1_constr = FALSE, family = "gaussian", seed = 2016)

write.csv(results$final_weights, file='weights.csv')