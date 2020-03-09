####################################################################
###############  Joint analysis of all plates    ###################
####################################################################


#### ---- Combine results into one matrix ----

# set location
topfolderlocation <- "xxxxxxxxxxxxxxx"

siRNApos <- read.table("TL_siRNA-Positions.csv", 
                       sep = ";", nrows = 40, skip = 4)

SummaryMatrix <- as.data.frame(read.csv("TL00931_AnalysisOutput/SummaryMatrix.csv", row.names = 1))
FeatureName <- colnames(SummaryMatrix)[1:23]
UpLim <- c(Inf,Inf,Inf,0.4,1e5,1e5,1e5,1e5,1e5,2e6,
           1e5,1e4,Inf,1.1,Inf,Inf,Inf,Inf,Inf,Inf,
           Inf,Inf,2e6)
LoLim <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,
           -Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,
           -Inf,-Inf,-Inf)
Limits <- data.frame(FeatureName = FeatureName,
                     UpLim = UpLim,
                     LoLim = LoLim)

PlateIDs <- c("00931", "00941", "00951", "00961", "00971", "00981")
siRNAposColIDs <- c(3,6,9,12,15,18)   # column in siRNApos which relates to the plate, in same order as in PlateIDs


for (i in 1:length(PlateIDs)) {
  filepath <- sprintf("TL%s_AnalysisOutput", PlateIDs[i])
  SummaryMatrices_add <- as.data.frame(read.csv(paste0(filepath, "/SummaryMatrix.csv"), row.names = 1))
  
  # Remove growth rate - associated data for 00941
  if (PlateIDs[i] == "00941" ) {
    SummaryMatrices_add$GrowthRate <- rep(NA, dim(SummaryMatrices_add)[1])
    SummaryMatrices_add$Compactness <- rep(NA, dim(SummaryMatrices_add)[1])                                   
    SummaryMatrices_add$Convexity <- rep(NA, dim(SummaryMatrices_add)[1])
    SummaryMatrices_add$Sphericity <- rep(NA, dim(SummaryMatrices_add)[1])
    SummaryMatrices_add$SpheroidGrowthRate <- rep(NA, dim(SummaryMatrices_add)[1])
    SummaryMatrices_add$SpheroidVolume <- rep(NA, dim(SummaryMatrices_add)[1])
    SummaryMatrices_add$CellVolume <- rep(NA, dim(SummaryMatrices_add)[1])
    
  }
  
  SummaryMatrices_add$PlateID <- PlateIDs[i]
  SummaryMatrices_add$Treatment <- siRNApos[2:39,siRNAposColIDs[i]]
  SummaryMatrices_add$siRNA <- siRNApos[2:39,siRNAposColIDs[i]+1]
  bin = matrix(NA, ncol=1, nrow=length(unique(SummaryMatrices_add$Treatment)))
  for (j in 1:length(unique(SummaryMatrices_add$Treatment))){
    bin[SummaryMatrices_add$Treatment==unique(SummaryMatrices_add$Treatment)[j]] <- 
      1 + (SummaryMatrices_add$siRNA[SummaryMatrices_add$Treatment==unique(SummaryMatrices_add$Treatment)[j]]==
             SummaryMatrices_add$siRNA[SummaryMatrices_add$Treatment==unique(SummaryMatrices_add$Treatment)[j]][1])
    
  }
  
  rownames(SummaryMatrices_add) <- paste0(PlateIDs[i], "_", rownames(SummaryMatrices_add), "_", SummaryMatrices_add$Treatment)
  
  if (i==1){
    FullMatrix <- SummaryMatrices_add
  }
  else {
    FullMatrix <- rbind(FullMatrix, SummaryMatrices_add)
  }
}

# Replace ESPLT = INCENP and add siRNA names for ESPL1 (s7423) and PLK1 (s448)

treats = as.character(FullMatrix$Treatment)
for (i in 1:dim(FullMatrix)[1]){
  treats[i] <- gsub("ESPL1", "INCENP", x=treats[i]) 
  rownames(FullMatrix) <- gsub("ESPL1", "INCENP", x=rownames(FullMatrix)) 
}

sirs <- as.character(FullMatrix$siRNA)
sirs[treats == "INCENP"] <- "s7423"
sirs[treats == "PLK1"] <- "s448"

FullMatrix$Treatment <- as.factor(treats)
FullMatrix$siRNA <- as.factor(sirs)

save(FullMatrix, file=paste0(filepaths.out, "/FullMatrix_no00941.Rdata"))


#### ---- Exclude outliers ----

FullMatrix_exclOutl <- FullMatrix[,1:length(UpLim)]
for (i in 1:length(UpLim)){
  isnotNA <- setdiff(1:dim(FullMatrix_exclOutl)[1], which(is.na(FullMatrix_exclOutl[,i])))
  isOutl <- which(FullMatrix_exclOutl[isnotNA,i]>UpLim[i])
  FullMatrix_exclOutl[isnotNA[isOutl],i] <- NA    # 39 outlier values excluded
}


#### ---- Heatmaps ---- 

# average all positions corresponding to each target and exclude empty positions

treats = as.character(unique(FullMatrix$Treatment))

FullMatrix_TargetAverage <- matrix(NA, ncol=ncol(FullMatrix_exclOutl), nrow = length(treats))
rownames(FullMatrix_TargetAverage) <- treats
colnames(FullMatrix_TargetAverage) <- colnames(FullMatrix_exclOutl)

for (i in 1:length(treats)){
  FullMatrix_TargetAverage[i,] <- colMeans(FullMatrix_exclOutl[FullMatrix$Treatment==treats[i],], na.rm=T)
}

FullMatrix_TargetAverage <- FullMatrix_TargetAverage[!(rownames(FullMatrix_TargetAverage)=="empty"),]

# Heatmaps

filepaths.out <- paste0(topfolderlocation, "/CombinedPlates")
hm.col <- colorRamps::matlab.like2(101)

# Rank-based clustering of summary matrix 
RankMatrix_TargetAverage <- FullMatrix_TargetAverage[,1:23]  
for (i in 1:dim(RankMatrix)[2]){   # rank entries: 
  RankMatrix_TargetAverage[,i] <- rank(RankMatrix_TargetAverage[,i], na.last = "keep", ties.method = "average")
}

# (i) based on Euclidean distance
d_TargetAverage <- dist(RankMatrix_TargetAverage, method = "euclidean") # distance matrix
fit_TargetAverage <- hclust(d_TargetAverage, method="complete") 
plot(fit_TargetAverage) # display dendogram
groups_TargetAverage <- cutree(fit_TargetAverage, k=5) # cut tree into 5 clusters
# draw dendogram with blue borders around the 5 clusters 
rect.hclust(fit_TargetAverage, k=5, border="blue")

heatmap.2(log(as.matrix(d_TargetAverage)/10000+1), trace="none")

# (ii) based on correlation
dd_TargetAverage <- as.dist((1 - cor(t(RankMatrix_TargetAverage)))/2)
dd_TargetAverage <- as.dist(sqrt(2*(1 - cor(t(RankMatrix_TargetAverage)))))

round(1000 * dd_TargetAverage) # (prints more nicely)
plot(hclust(dd_TargetAverage)) # to see a dendrogram of clustered variables

fit_TargetAverage <- hclust(dd_TargetAverage, method="complete") 
plot(fit_TargetAverage) # display dendogram
groups_TargetAverage <- cutree(fit_TargetAverage, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit_TargetAverage, k=6, border="green")

heatmap.2(log(as.matrix(dd_TargetAverage)/10000+1), trace="none")

# save heatmap

pdf(file = paste0(filepaths.out, "/Heatmap_complete_scaledCol.pdf"),
    width = 11.8, height = 14.2, pointsize = 18)
aheatmap(RankMatrix_TargetAverage,
         distfun = "euclidean", hclustfun = "complete",
         annCol = groups_TargetAverage,
         scale = "col")
dev.off()

