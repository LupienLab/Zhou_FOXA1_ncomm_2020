setwd("/Users/lupienlab/Desktop/lupien_lab_stanley/PCa_Project/AchillesData")

### This reads in the achilles data
achilles<-read.table("Achilles_QC_v2.4.3.rnai.Gs.gct", skip=2, header=TRUE, sep="\t")
ach<-achilles[,c(2,3,216)]

### Get your genes of interest
GenesOfInterest<-read.table("/Users/lupienlab/Desktop/lupien_lab_stanley/PCa_Project/AchillesData/GenesOfInterest.txt", sep="\t", stringsAsFactors = FALSE)[,1]

### Find the number of genes of interest in the Achilles data
wach<-c()
for(i in 1:length(GenesOfInterest)) {
  wach<-c(wach, grep(pattern=GenesOfInterest[i], x=ach[,1]))
}

#wach<-wach[-c(1,3,4,6,17,18,19,23,24,25,26,27,28)]

floop<-1:nrow(ach)[-wach]
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
  pMat[i,1]<-mean(ach[sample(floop, SizeOfDataset),2])
  pMat[i,2]<-mean(ach[sample(floop, SizeOfDataset),3])
}

Enri<-c(mean(ach[wach,2]), mean(ach[wach,3]))

TestEnrichment(pMat, Enri, c("22Rv1","VCaP"), c("#fdae61","#2c7bb6"), "Achilles Data", 1.959964)

geom_text(data=filter(results, padj<0.05), aes(label=Gene))