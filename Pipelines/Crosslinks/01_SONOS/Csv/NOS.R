

setwd('C:/Dev/Github/LeucipPipelines/Pipelines/Crosslinks/01_SONOS/Csv/')
ddatanos <- read.csv("NOS_MEDIUM_03_Geometry.csv",stringsAsFactors = FALSE, row.names=NULL)

ddatanos$NZ.SG = as.numeric(ddatanos$NZ.SG) 
ddatanos$NZ.O = as.numeric(ddatanos$NZ.O) 
ddatanos$NZ.N = as.numeric(ddatanos$NZ.N) 
ddatanos$NZ.SG2 = as.numeric(ddatanos$NZ.SG2) 
ddatanos$NZ.NO2 = as.numeric(ddatanos$NZ.NO2) 
ddatanos$NZ.NO3 = as.numeric(ddatanos$NZ.NO3) 

numericdatanos<- dplyr::select (ddatanos,NZ.SG,NZ.O,NZ.N,NZ.SG2,NZ.NO2,NZ.NO3)

matrix_datanos <- cbind(ddatanos$NZ.SG,ddatanos$NZ.O,ddatanos$NZ.N,ddatanos$NZ.SG2,ddatanos$NZ.NO2,ddatanos$NZ.NO3) 

len <- nrow(ddatanos)

grouped_datanos <- data.frame(matrix_datanos,as.factor(c(rep("A",len),rep("B",len),rep("C",len),rep("D",len),rep("E",len),rep("F",len)))) 


rows = nrow(matrix_datanos)

pcanos <- prcomp(matrix_datanos, center = TRUE,scale=FALSE) 
plot(pcanos$x[, 1], pcanos$x[, 2], pch=21, bg=c("blue","cyan","green","yellow","black","pink")[grouped_datanos[,7]], main = "PCA", xlab = "PC1", ylab = "PC2") 

#K means analysis

mydatanos<-data.frame(pcanos$x[, 1],pcanos$x[, 2],pcanos$x[, 3],pcanos$x[, 4])  

withingroupssnos <- (nrow(mydatanos)-1)*sum(apply(mydatanos,2,var)) 
for (i in 2:8) withingroupssnos[i] <- sum(kmeans(mydatanos,centers=i)$withinss) 
par(pin=c(2,2),font=2,ps=10,family="sans") 
plot(1:8, withingroupssnos[1:8], type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") 

# 3 looks like a good number
myclustersnos <- kmeans(mydatanos,3) 

library(cluster) 
clusplot(mydatanos, myclustersnos$cluster,color=TRUE, shade=TRUE, labels=1, lines=0) 
library(cluster) 

#(need to press ESC with mouse in graphics window to get out of clusplot) 

