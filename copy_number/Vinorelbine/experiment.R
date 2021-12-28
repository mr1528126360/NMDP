source('CLSPCA.R')
library(glmnet)
myString <- "read seq data!"
print(myString)
data <- read.csv("cnv_zen.csv", header=TRUE)
myString <- "read pathway data!"
print(myString)
RelationMatrix <- read.csv("../pathway.csv", header=FALSE)
myString <- "start!"
print(myString)
edges = list()
gene1 = data.matrix(data)
gene = t(gene1)
for (i in 1:317757){
  edges[[i]] = c((RelationMatrix[i,1]+1),(RelationMatrix[i,2]+1))
}
n = 5000
out3 =  CLSPCA(gene, k = 3, edges, k.group=n,we=0.3, t = 0.1, niter=10, err=0.0001, Num.init=5)
write.csv (out3[["V"]], file ="CLSPCA_V_CNV.csv")
write.csv (out3[["U"]], file ="CLSPCA_U.csv")
myString <- "end!"
print(myString)
