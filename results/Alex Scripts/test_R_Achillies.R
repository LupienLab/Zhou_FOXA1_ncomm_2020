setwd("/users/lupienlab2/Desktop/Prostate_cancer_Achillies")
AchilliesCRISPR <-read.table("Achilles_v3.3.8.Gs.gct", skip=2, header=TRUE, sep="\t")
interest <-AchilliesCRISPR [,c(2,18,28)]

GenesOfInterest<-read.table("/users/lupienlab2/Desktop/Prostate_cancer_Achillies/Core_PRC1.txt", sep="\t", stringsAsFactors = FALSE)[,1]

wach<-c()
for(i in 1:length(GenesOfInterest)) {
  wach<-c(wach, grep(pattern=GenesOfInterest[i], x=interest[,1]))
}
#wach <- wach[-c(4,6,31,8,9,32,33,34,35,43,44,45,46,47,48,49,50,51,52,53,54,55,56)]

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
  for(i in 1:ncol(MedPerm)) {MedPerm[,i]<-(Permutations[,i]-Means[i])/SDs[i]}
  for(i in 1:ncol(MedPerm)) {EnNorm[i]<-(EnNorm[i]-Means[i])/SDs[i]}
  boxplot(MedPerm, las=2, ylab="Mean Normalized Achilles Gene Level Score", cex.axis=0.8, names=CellLines, 
          col=DesiredColours, main=title, ylim=c(-6,6))
  csc<-rep("red", length(EnNorm))
  csc[which(EnNorm<0)]<-"blue"
  csc[which(abs(EnNorm)<Corrected)]<-"black"
  points(1:length(EnNorm), EnNorm, pch=18, cex=4, col="red")  
  abline(h=-Corrected, col="red", lty=2)
  abline(h=Corrected, col="red", lty=2)
}


DesiredNumberOfPermutations<-1000
pMat<-matrix(ncol=2, nrow=DesiredNumberOfPermutations)
for(i in 1:DesiredNumberOfPermutations) {
  pMat[i,1]<-mean(interest[sample(floop, SizeOfDataset),2])
  pMat[i,2]<-mean(interest[sample(floop, SizeOfDataset),3])
}

Enri<-c(mean(interest[wach,2]), mean(interest[wach,3]))

TestEnrichment(pMat, Enri, c("LNCAP","PC-3"), c("#fdae61","#2c7bb6"), "Achilles Data", 1.959964)

