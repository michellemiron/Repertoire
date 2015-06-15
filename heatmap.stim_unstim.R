#Make heatplot of alloreactive clones
#Plot Change between Stim and Unstim samples for CD8 OR CD4
  #Example inputfile
    #aminoAcid	vGeneName	jGeneName	total	stim_CD4	stim_CD8	unstim_CD4	unstim_CD8
    #CASIGPLNEKLFF	TCRBV25-01	TCRBJ01-04	243169	619	204005	33	38512

#Run in R by:
      # heatmap("filename",column= cd4change, titleplot ="cd4change")

##Color Key:
#heatmap of CD4freq in stim - CD4freq in unstim / (CD4 freq in stim + CD4 freq in unstim)
  #grey : not present in either
  #red = 1 : only present in stim
  #blue = 1 : only present i unstim
  #red > 0 : higher in stim  
  #blue <0 : higher in unstim 
      #alloreactive is red >= 0.66 : 5 fold expansion = >= 0.66 


cleanup <- function(cd4, cd8, ratio = 5) {
  
  ## remove contaminated clones
  
  ambi = (cd4 > 0 & cd8 > 0 & cd4 / cd8 > 1/ratio & cd4 / cd8 < ratio) 
  cd4exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 <= 1/ratio ) 
  cd8exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 >= ratio )
  #throws away (marks as true for exlcusion) ones that are ambi in stim, but not ambi in unstim
  #if you want to keep cd8 that is ambi in stim with cd4, but clearly not ambi in unstim and dominant, need to do comparison against all
  #remember, true means to be excluded and false means keep. false could be a zero so length of data==F is not number to keep
  # print(paste(length(ambi[ambi==T]), length(cd4exclude[cd4exclude==T]), length(cd8exclude[cd8exclude==T])))
  
  return(cbind(ambi | cd4exclude, ambi | cd8exclude))
  
}

threshold <- function(stim, unstim, ratio = 5) {
  
  ## remove clones under the threshold
  
  alloeexclude = (stim > 0 & unstim > 0 & stim / unstim <= 5 ) 
  
  return(alloexclude)
  
}

normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i] / sum(data[,i])
  }
  return(data)
}

##plot VJ by aggregate data

heatmap <- function(file, clean=T, threshold=F, column= cd8change, titleplot ="cd8change") {
  file <- "HC_19_10_count_aa.tsv"

setwd("/Volumes/BigData/TCR/Sykes/alloresponse_count_files copy/")
#get and read file
a = read.table(file, header=T)

### Set Data
cd4reads = cbind(a[,5], a[,7])
cd8reads = cbind(a[,6], a[,8])

#(stim, unstim)  
#normalize function
cd4 = normalize(cd4reads)
cd8 = normalize(cd8reads)

cleans  = "raw"

if (clean == T) {
#True False table of overlapping clones. True is contimanation to be removed
rows1 = cleanup(cd4[,1], cd8[,1]) #stim
rows2 = cleanup(cd4[,2], cd8[,2]) #unstim
m1 = cbind(rows1,rows2)
m1 = 1-m1 #TRUE = 0 and FALSE = 1
m2 = cbind(cd4[,1], cd8[,1],cd4[,2], cd8[,2]) #reference data
m3 = m1*m2 
colnames(m3) <- c("cd4stim","cd8stim","cd4unstim","cd8unstim")
#m3 is clean data with zeros placed where overlapping clones were

cd4 = cbind(m3[,1],m3[,3])
cd8 = cbind(m3[,2], m3[,4])                                                                                                                                                                                                                                            
cd4 = normalize(cd4)
cd8 = normalize(cd8)
}

cleans = "clean"

##aggregate by VJ genes

v <- a[,2]
j <- a[,3] 
cd4agg <- aggregate(x=cd4, by= list(j,v), FUN= "sum")
cd8agg <- aggregate(x=cd8, by= list(j,v), FUN= "sum")
cd4 <- cd4agg[,3:4]
cd8 <- cd8agg[,3:4]
colnames(cd4agg) <- c("X1", "X2", "value", "value")
colnames(cd8agg) <- c("X1", "X2", "value", "value")


cd4stim <- cd4agg[,1:3]
cd4unstim <- cd4agg[,-3]
cd4change <- (cd4stim[,3]-cd4unstim[,3])/(cd4stim[,3]+cd4unstim[,3])
cd4change<- cbind(cd4stim[,-3], cd4change)
colnames(cd4change) <- c("X1", "X2","value")
cd8stim <- cd8agg[,1:3]
cd8unstim <- cd8agg[,-3]
cd8change <- (cd8stim[,3]-cd8unstim[,3])/(cd8stim[,3]+cd8unstim[,3])
cd8change<- cbind(cd8stim[,-3], cd8change)
colnames(cd8change) <- c("X1", "X2","value")


column = column
title <- paste( titleplot, file, sep= " ")


##simple plot of frequency
library(ggplot2)
p <- ggplot(column, aes(y=X2, x=X1))
p + geom_tile(aes(fill=value)) + scale_fill_gradient2(mid="black", high="red", low="blue") + xlab("") + ylab("") + ggtitle(title)

}
#scale_fill_gradient(low="blue", high="red") 


