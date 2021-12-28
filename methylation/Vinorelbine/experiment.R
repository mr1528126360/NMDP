# -*- coding: utf-8 -*-
source('CLSPCA.R')
library(glmnet)
myString <- "read seq data!"
print(myString)
#名称要对应
data <- read.csv("Methylation_zen.csv", header=TRUE)
myString <- "read pathway data!"
print(myString)
RelationMatrix <- read.csv("../pathway.csv", header=FALSE)
myString <- "start!"
print(myString)
edges = list()
gene1 = data.matrix(data)
gene = t(gene1)
for (i in 1:139017209){
  edges[[i]] = c((RelationMatrix[i,1]+1),(RelationMatrix[i,2]+1))
}
n = 5000
out3 =  CLSPCA(gene, k = 2, edges, k.group=n,we=0.3, t = 0.1, niter=10, err=0.0001, Num.init=3)
write.csv (out3[["V"]], file ="CLSPCA_V_Methylation.csv")
write.csv (out3[["U"]], file ="CLSPCA_U_Methylation.csv")
myString <- "end!"
print(myString)
