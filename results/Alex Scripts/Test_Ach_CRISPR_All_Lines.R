library("ggplot2")

setwd("/Users/lupienlab/Desktop/lupien_lab_stanley/PCa_Project/AchillesData")
AchilliesCRISPR <-read.table("Achilles_v3.3.8.Gs.gct", skip=2, header=TRUE, sep="\t")
interest <-AchilliesCRISPR [,c(2,18,28)]
interest <-AchilliesCRISPR [,-c(1)]
GenesOfInterest<-read.table("/Users/lupienlab/Desktop/lupien_lab_stanley/PCa_Project/AchillesData/GenesofInterest.txt", sep="\t", stringsAsFactors = FALSE)[,1]

wach<-c()
for(i in 1:length(GenesOfInterest)) {
  wach<-c(wach, grep(pattern=GenesOfInterest[i], x=interest[,1]))
}
wach<-wach[-c(4,6,31,32,33,34,35,43,44,45,46,47,48,49,50,51,52,53,54,55,56)] 


floop<-1:nrow(interest)[-wach]
SizeOfDataset<-length(wach)

sloop<-sample(floop, SizeOfDataset)

TestEnrichment<-function(Permutations, Enrichments, CellLines, DesiredColours, title, Corrected) {
  #Permutations = Matrix where permutations are in rows, distinct entities being tested are in columns
  #Enrichments = Vector of ncol(Permutations) where the observed value is shown
  #CellLines = Names of parameters
  #DesiredColours = Colour of boxplots
  #title = Title of Graph
  #Corrected = Z-Score for significance.
  
  Means<-apply(Permutations, 2, mean)
  SDs<-apply(Permutations, 2, sd)
  MedPerm<-Permutations
  EnNorm<-Enrichments
  # normalize Permutations
  for(i in 1:ncol(MedPerm)) {MedPerm[,i]<-(Permutations[,i]-Means[i])/SDs[i]}
  # normalize Enrichments
  for(i in 1:ncol(MedPerm)) {EnNorm[i]<-(EnNorm[i]-Means[i])/SDs[i]}
  
  # generate boxplot
  boxplot(
    MedPerm,
    las=2,
    ylab="Mean Normalized Achilles Gene Level Score",
    cex.axis=0.8,
    names=CellLines, 
    col=DesiredColours,
    main=title,
    ylim=c(-6,6),
    cex.names=0.5
  )
  
  # colour dots of interest
  # red by default
  csc<-rep("red", length(EnNorm))
  # blue if "essential"
  csc[which(EnNorm<0)]<-"blue"
  # black if not significant
  csc[which(abs(EnNorm)<Corrected)]<-"black"
  # generate points on top of boxplot
  points(
    1:length(EnNorm),
    EnNorm,
    pch=18,
    cex=4,
    col="red"
  )
  abline(h = Corrected, col = "red", lty = 2)
  abline(h = -Corrected, col = "red", lty = 2)
 
}


DesiredNumberOfPermutations<-1000
pMat<-matrix(ncol=2, nrow=DesiredNumberOfPermutations)
for(i in 1:DesiredNumberOfPermutations) {
  pMat[i,1]<-mean(interest[sample(floop, SizeOfDataset),2])
  pMat[i,2]<-mean(interest[sample(floop, SizeOfDataset),3])
}

# interest[, 17] is LNCaP, interest[, 27] is PC3
Enri<-c(interest[wach, 17], interest[wach, 27])
colnames(pMat) <- colnames(Enri)

# multiple test corrected z-score that would classified as significant
sig_zscore <- 1.959964
TestEnrichment(pMat, Enri, c("LNCaP","PC3"), c("#fdae61","#2c7bb6"), "Achilles Data", sig_zscore)


DesiredNumberOfPermutations<-10000
pMat<-matrix(ncol=ncol(interest)-1, nrow=DesiredNumberOfPermutations)
for(i in 1:DesiredNumberOfPermutations) {
  for(j in 1:ncol(pMat)) {
    pMat[i,j]<-mean(interest[sample(floop, SizeOfDataset),j+1])
  }
}

Enri<-apply(data.matrix(interest[wach,-c(1)]), MARGIN = 2, mean)

TestEnrichment(pMat[,order(Enri)], Enri[order(Enri)], colnames(interest)[2:ncol(interest)][order(Enri)], 
               c("red","blue","darkorchid3","darkgrey","blue","darkorchid3","darkgrey","yellow","red","green","darkgrey","darkorange","darkgrey","saddlebrown","pink1","darkorange","darkgrey","darkorange",
                 "darkgrey","yellow","darkgrey","darkorchid3","darkorange","darkorange","lightseagreen","yellow","darkorchid3","darkorange","darkgrey","darkorange","darkorange", "darkgrey", "yellow"),
              # rep("lightsalmon1", ncol(interest)-1), 
               "Achilles Data", sig_zscore)