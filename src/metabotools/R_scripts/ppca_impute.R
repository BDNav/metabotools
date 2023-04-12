

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
message(args)   

library("MetaboAnalystR")
mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, args[1], "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeMissingVar(mSet, method="ppca")
mSet<-SaveTransformedData(mSet)