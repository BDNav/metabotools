
# requires NMF library:
# install.packages('NMF')


library(data.table)
library(NMF)

args <- commandArgs(trailingOnly = TRUE)

data <- fread(args[1])
data2<-data[,-c(1:2)]
data3<-t(data2)

res <- nmf(data3, 4) # 4 is hardcoded for now

fit(res)


V.hat <- fitted(res)
dim(V.hat)
 
summary(res)


# get matrix W
w <- basis(res)
dim(w)
# get matrix H
h <- coef(res)
dim(h)


s <- featureScore(res)
summary(s)
# compute the scores and characterize each metagene
s <- extractFeatures(res)
str(s)

write.csv(h, file='coefs.csv')
write.csv(w, file='basis.csv')
# fileConn<-file("s.txt")
# writeLines(str(s), fileConn)
# close(fileConn)

# cat(str(s), file='s.txt')

writeLines(sapply(s, toString), "s.txt")

#res.multirun <- nmf(esGolub, 3, nrun=5)



#estim.r <- nmf(esGolub, 2:6, nrun=10, seed=123456)