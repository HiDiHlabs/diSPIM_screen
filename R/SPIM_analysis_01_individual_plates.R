### Bjoern data analysis - final ###

#### ---- Load libraries ----

library(ggplot2)
library(data.table)
library(reshape2)
library(gplots)
library(grid)
library(colorRamps)
library(NMF)

# set location
topfolderlocation <- "xxxxxxxxxxxxxxx"

siRNApos <- read.table(sprintf("%s/TL_siRNA-Positions.txt",topfolderlocation),
                       sep = "\t", nrows = 40, skip = 4)

#PlateIDs <- c("00931", "00941", "00951", "00961", "00971", "00981")
#siRNAposColIDs <- c(3,6,9,12,15,18)   # column in siRNApos which relates to the plate, in same order as in PlateIDs
PlateIDs <- "00971"
siRNAposColIDs <- 15


for (i in 1:length(PlateIDs)) {  # 1:length(PlateIDs)
  
  filepaths <- sprintf("%s/TL%s_S_hSPIM_Output", topfolderlocation, PlateIDs[i])
  filepaths.out <- sprintf("%s/TL%s_AnalysisOutput", topfolderlocation, PlateIDs[i])
  
  ############################################################
  ########## Analysis of spheroid properties        ##########
  ############################################################
  
  SpheroidPropertiesIn <- list()
  for (pos in 1:38){
    SpheroidPropertiesIn[[pos]] <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_SpheroidProperties/Pos%s.csv", topfolderlocation, PlateIDs[i], pos),
                                              sep = ",", nrows = 290, skip = 1)
    
  }
  
  SpheroidProperties <- list()
  for (pos in 1:38){
    SpheroidProperties[[pos]] = cbind(SpheroidPropertiesIn[[pos]][,c(9,3,4,5,6)])
    colnames(SpheroidProperties[[pos]]) <- c("Time Point", "Spheroid_Volume", "Spheroid_Compactness", 
                                             "Spheroid_Convexity", "Spheroid_Sphericity")
  }
  
  
  ## Compactness / Convexity / Sphericity
  ############################################################
  
  FeatureName = c("Spheroid_Compactness", "Spheroid_Convexity", "Spheroid_Sphericity")
  DeathTimepoints = matrix(nrow = length(FeatureName), ncol = 38)
  rownames(DeathTimepoints) <- FeatureName
  
  for (f in 1:3){
    
    FeatureIdx <- f
    
    # user-defined function to label data frame columns by position
    posfunc <- function(x){
      y <- sprintf("Pos_%s", x)
      return(y)
    }
    
    SpheroidFeat <- matrix(NA, ncol=38, nrow=290)  # initialise
    colnames(SpheroidFeat) <- sapply(X = seq(1,38,1), FUN = posfunc)
    
    for (pos in 1:38){
      SpheroidFeat[,pos] <- SpheroidProperties[[pos]][,FeatureName[FeatureIdx]]
      if(sum(is.na(SpheroidFeat[,pos]))>0){      # cell is dead
        DeathTimepoints[FeatureIdx,pos] <- min(which(is.na(SpheroidFeat[,pos])))
      }
    }
    
    # plot
    
    SpheroidFeat_toplot <- data.frame(SpheroidFeat)
    SpheroidFeat_toplot$TimePoint <- seq(1,290,1)
    
    pdf(sprintf(paste0(filepaths.out, "/%s", "_TimeSeries.pdf"), FeatureName[FeatureIdx]),width=15,height=4,paper='special') 
    print(ggplot(data=melt(SpheroidFeat_toplot, id="TimePoint"),
                 aes(x=TimePoint, y=value, colour=variable)) +
            geom_line())
    dev.off()
    
    # violin plot across all time points
    pdf(sprintf(paste0(filepaths.out, "/%s", "_VlnPlot.pdf"), FeatureName[FeatureIdx]),width=15,height=4,paper='special') 
    print(ggplot(data=melt(SpheroidFeat_toplot, id="TimePoint"),
                 aes(x=variable, y=value, fill=variable)) +
            geom_violin() +
            guides(fill=FALSE))
    dev.off()
    
  }
  
  
  ## Spheroid volume (based on hull)
  ############################################################
  
  SpheroidVol <- matrix(nrow = 290, ncol = 38)
  # user-defined function to label data frame columns by position
  posfunc <- function(x){
    y <- sprintf("Pos_%s", x)
    return(y)
  }
  colnames(SpheroidVol) <- sapply(X = seq(1,38,1), FUN = posfunc)
  SpheroidVol <- data.frame(SpheroidVol)
  
  for (pos in 1:38){
    SpheroidVol[,pos] <- SpheroidProperties[[pos]][,"Spheroid_Volume"]
  }
  
  # plot
  
  SpheroidVol_toplot <- data.frame(SpheroidVol)
  SpheroidVol_toplot$TimePoint <- seq(1,290,1)
  # replace NaN with 0
  SpheroidVol_toplot[is.na(SpheroidVol_toplot)] <- 0
  
  
  pdf(paste0(filepaths.out, "/SpheroidVolume_TimeSeries.pdf"),width=15,height=4,paper='special') 
  print(ggplot(data=melt(SpheroidVol_toplot, id="TimePoint"),
               aes(x=TimePoint, y=value, colour=variable)) +
          geom_line() +
          ylim(0, 2000000))
  dev.off()

  ############################################################
  ########## Analysis of nuclei properties          ##########
  ############################################################
  

  ### Bar charts of % of time spent in different phases (across all images)
  ############################################################
  
  ClassPercentage <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_NucleiProperties/ClassPercentage.csv", topfolderlocation, PlateIDs[i]),
                                sep = ",", nrows = 152, skip = 1,
                                col.names = c("RowID",	"CellPhase",	"Pos",	"PosNr",	"Rank",	"FractionPhase"))
  
  df <- ClassPercentage
  
  PhaseNames <- c("Interphase", "Prophase", "Metaphase", "Anaphase")
  PhaseNames_initials <- c("I", "P", "M", "A") 
  
  ClassPercentageBarplots <- list()
  for (j in 1:38){
    dfsub <- subset(df,PosNr==j)
    p1 <- ggplot(data = dfsub[match(PhaseNames, dfsub$CellPhase),], aes(x=CellPhase, y=FractionPhase, fill=CellPhase)) +
      geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_brewer(palette="Set3") +
      theme(legend.position = "none") +
      scale_x_discrete(breaks=PhaseNames,labels=PhaseNames_initials) +
      ylim(0, 100)
    ClassPercentageBarplots[[j]] <- p1  # add each plot into plot list
  }
  
  pdf(paste0(filepaths.out, "/ClassPercentageBarplots.pdf"),width=15,height=12,paper='special') 
  print(multiplot(plotlist = ClassPercentageBarplots, cols = 8)) # <-- slow/large
  dev.off()
  
  # Rearrange data frame
  df <- reshape(ClassPercentage[,c("PosNr","CellPhase","FractionPhase")],timevar="CellPhase",idvar="PosNr",direction="wide")
  pdf(paste0(filepaths.out, "/FractionInInterphase_Histogram.pdf"),width=5,height=5,paper='special') 
  print(hist(df$FractionPhase.Interphase))
  dev.off()
  
  
  ### Number of nuclei and growth rate
  ############################################################
  
  NucleiPropertiesIn <- list()
  for (pos in 1:38){
       NucleiPropertiesIn[[pos]] <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_NucleiProperties/Pos%s.csv", topfolderlocation, PlateIDs[i], pos),
                                            sep = ",", nrows = 290, skip=1)
  }
  
  NucleiProperties <- list()
  for (pos in 1:38){
    NucleiProperties[[pos]] = cbind(NucleiPropertiesIn[[pos]][,c(3,5,6,7,8,9,13,14,15,16,29,34,39,44)])
    colnames(NucleiProperties[[pos]]) <- c("TimePoint", "Total_Count", 
                                           "Interphase_Count",	"Prophase_Count",	"Metaphase_Count",	"Anaphase_Count",
                                           "Interphase_MedianVoxels", "Prophase_MedianVoxels", "Metaphase_MedianVoxels", "Anaphase_MedianVoxels",
                                           "Interphase_STDVoxels", "Prophase_STDVoxels", "Metaphase_STDVoxels", "Anaphase_STDVoxels"
    )
  }
  
  # not all positions have the same number of time points -> determine n.rows for each
  n.rows <- list()
  for (pos in 1:38){
    n.rows[[pos]] <- dim(NucleiProperties[[pos]])[1]
  }
  
  
  ## Number of cells
  
  NrCells <- matrix(NA, nrow = 290, ncol = 38)
  # user-defined function to label data frame columns by position
  posfunc <- function(x){
    y <- sprintf("Pos_%s", x)
    return(y)
  }
  colnames(NrCells) <- sapply(X = seq(1,38,1), FUN = posfunc)
  NrCells <- data.frame(NrCells)
  
  for (pos in 1:38){
    NrCells[1:dim(NucleiProperties[[pos]])[1], pos] <- NucleiProperties[[pos]][,"Total_Count"]
  }
  
  # plot
  
  NrCells_toplot <- data.frame(NrCells)
  NrCells_toplot$TimePoint <- seq(1,290,1)
  
  pdf(paste0(filepaths.out, "/NrNuclei_TimeSeries.pdf"),width=15,height=12,paper='special') 
  ggplot(data=melt(NrCells_toplot, id="TimePoint"),
         aes(x=TimePoint, y=value, colour=variable)) +
    geom_line()
  dev.off()
  
  
  ## Growth rate (based on smoothed data)
  GrowthRate <- list()
  for (pos in 1:38){
    GrowthRate[[pos]] <- data.frame(Nr_LOWESS=numeric(length = n.rows[[pos]]),
                                    InstantGrowthRate=numeric(length = n.rows[[pos]]),
                                    stringsAsFactors=FALSE)
    GrowthRate[[pos]]$Nr_LOWESS <- lowess(NrCells[1:n.rows[[pos]],pos], f=1/3, iter=3L, delta=0.01 * diff(range(NrCells[1:n.rows[[pos]],pos])))$y
    GrowthRate[[pos]]$InstantGrowthRate[2:n.rows[[pos]]] <- diff(GrowthRate[[pos]]$Nr_LOWESS)
  }
  
  
  # plot
  
  GrowthRates_toplot <- matrix(nrow=289,ncol=38)
  colnames(GrowthRates_toplot) <- sapply(X = seq(1,38,1), FUN = posfunc)
  for (pos in 1:38){
    GrowthRates_toplot[,pos] <- GrowthRate[[pos]]$InstantGrowthRate[2:290]
  }
  GrowthRates_toplot <- data.frame(GrowthRates_toplot)
  GrowthRates_toplot$TimePoint <- seq(1,289,1)
  
  pdf(paste0(filepaths.out, "/GrowthRates_TimeSeries.pdf"),width=15,height=12,paper='special') 
  ggplot(data=melt(GrowthRates_toplot, id="TimePoint"),
         aes(x=TimePoint, y=value, colour=variable)) +
    geom_line()
  dev.off()
  
  
  ## LOWESS of instantaneous growth rates
  GrowthRate_smoothed <- list()
  for (pos in 1:38){
    GrowthRate_smoothed[[pos]] <- data.frame(
      InstantGrowthRate_LOWESS=numeric(length = dim(NucleiProperties[[pos]])[1]),
      stringsAsFactors=FALSE)
    GrowthRate_smoothed[[pos]]$InstantGrowthRate_LOWESS <- lowess(GrowthRate[[pos]]$InstantGrowthRate, f=1/3, iter=3L, delta=0.01 * diff(range(GrowthRate[[pos]]$InstantGrowthRate)))$y
  }
  
  # plot
  
  GrowthRate_smoothed_toplot <- matrix(NA,nrow=289,ncol=38)
  colnames(GrowthRate_smoothed_toplot) <- sapply(X = seq(1,38,1), FUN = posfunc)
  for (pos in 1:38){
    GrowthRate_smoothed_toplot[1:(n.rows[[pos]]-1),pos] <- GrowthRate_smoothed[[pos]]$InstantGrowthRate_LOWESS[2:n.rows[[pos]]]
  }
  GrowthRate_smoothed_toplot <- data.frame(GrowthRate_smoothed_toplot)
  GrowthRate_smoothed_toplot$TimePoint <- seq(1,289,1)
  
  pdf(paste0(filepaths.out, "/GrowthRates_TimeSeries_LOWESS.pdf"),width=15,height=12,paper='special') 
  ggplot(data=melt(GrowthRate_smoothed_toplot, id="TimePoint"),
         aes(x=TimePoint, y=value, colour=variable)) +
    geom_line()
  dev.off()
  
  
  ## Mean growth rate (based on LOWESS and linear fit)
  MeanGrowthRate <- matrix(nrow=38, ncol=4)
  colnames(MeanGrowthRate) <- c("Mean", "STD", "LinearFit_slope", "LinearFit_intercept")
  for (pos in 1:38){
    MeanGrowthRate[pos,1] <- mean(GrowthRate[[pos]]$InstantGrowthRate[-1])
    MeanGrowthRate[pos,2] <- sd(GrowthRate[[pos]]$InstantGrowthRate[-1])
    linfit <- lm(NrCells[1:n.rows[[pos]],pos] ~ seq(1,n.rows[[pos]],1))
    MeanGrowthRate[pos,3] <- linfit[[1]][2]
    MeanGrowthRate[pos,4] <- linfit[[1]][1]
  }
  
  # plot
  MeanGrowthRate_toplot <- data.frame(MeanGrowthRate)
  
  pdf(paste0(filepaths.out, "/MeanGrowthRates_Histogram.pdf"),width=5,height=4,paper='special') 
  ggplot(data=MeanGrowthRate_toplot,
         aes(x=Mean)) +
    geom_histogram(colour="black", fill="grey", position="dodge", binwidth = 0.02)
  dev.off()
  
  
  ### Size of nuclei in different phases (averaged)
  ############################################################
  
  VoxelsPerCell <- matrix(nrow=38, ncol=10)
  colnames(VoxelsPerCell) <- c("Interphase_mean", "Interphase_sd", "Prophase_mean", "Prophase_sd", "Metaphase_mean", "Metaphase_sd", "Anaphase_mean", "Anaphase_sd", "All_mean", "All_sd")
  
  for (pos in 1:38){
    VoxelsPerCell[pos,1] <- mean(NucleiProperties[[pos]]$Interphase_MedianVoxels[NucleiProperties[[pos]]$Interphase_MedianVoxels>0.1])
    VoxelsPerCell[pos,2] <- sd(NucleiProperties[[pos]]$Interphase_MedianVoxels[NucleiProperties[[pos]]$Interphase_MedianVoxels>0.1])
    VoxelsPerCell[pos,3] <- mean(NucleiProperties[[pos]]$Prophase_MedianVoxels[NucleiProperties[[pos]]$Prophase_MedianVoxels>0.1])
    VoxelsPerCell[pos,4] <- sd(NucleiProperties[[pos]]$Prophase_MedianVoxels[NucleiProperties[[pos]]$Prophase_MedianVoxels>0.1])
    VoxelsPerCell[pos,5] <- mean(NucleiProperties[[pos]]$Metaphase_MedianVoxels[NucleiProperties[[pos]]$Metaphase_MedianVoxels>0.1])
    VoxelsPerCell[pos,6] <- sd(NucleiProperties[[pos]]$Metaphase_MedianVoxels[NucleiProperties[[pos]]$Metaphase_MedianVoxels>0.1])
    VoxelsPerCell[pos,7] <- mean(NucleiProperties[[pos]]$Anaphase_MedianVoxels[NucleiProperties[[pos]]$Anaphase_MedianVoxels>0.1])
    VoxelsPerCell[pos,8] <- sd(NucleiProperties[[pos]]$Anaphase_MedianVoxels[NucleiProperties[[pos]]$Anaphase_MedianVoxels>0.1])
    a <- NucleiProperties[[pos]][,c("Interphase_MedianVoxels", "Prophase_MedianVoxels", "Metaphase_MedianVoxels", "Anaphase_MedianVoxels")]
    VoxelsPerCell[pos,9] <- mean(a[a>0.1])
    VoxelsPerCell[pos,10] <- sd(a[a>0.1])
  }
  
  # plot
  
  VoxelsPerCell_toplot <- data.frame(VoxelsPerCell)
  VoxelsPerCell_toplot$Pos <- seq(1,38,1)
  # replace NaN with 0
  VoxelsPerCell_toplot[is.na(VoxelsPerCell_toplot)] <- 0
  
  pdf(paste0(filepaths.out, "/VoxelsPerCell_VlnPlot.pdf"),width=10,height=5,paper='special') 
  ggplot(data=melt(VoxelsPerCell_toplot[,c(1,3,5,7,11)], id="Pos"),
         aes(x=variable, y=value, fill=variable)) +
    geom_violin(show.legend=F)
  dev.off()
  
  pdf(paste0(filepaths.out, "/VoxelsPerCell_Histogram.pdf"),width=10,height=5,paper='special') 
  ggplot(data=melt(VoxelsPerCell_toplot[,c(1,3,5,7,11)], id="Pos"),
         aes(x=value, colour=variable)) +
    geom_histogram(fill="white", position="dodge", binwidth = 50000)
  dev.off()
  
  
  ## Volume per nucleus (ratio spheroid volume to cell number)
  ############################################################
  
  CellVol <- matrix(NA, nrow = 290, ncol = 38)
  # user-defined function to label data frame columns by position
  posfunc <- function(x){
    y <- sprintf("Pos_%s", x)
    return(y)
  }
  colnames(CellVol) <- sapply(X = seq(1,38,1), FUN = posfunc)
  CellVol <- data.frame(CellVol)
  
  for (pos in 1:38){
    CellVol[1:n.rows[[pos]],pos] <- SpheroidProperties[[pos]][1:n.rows[[pos]],"Spheroid_Volume"]/
      NucleiProperties[[pos]][,"Total_Count"]
  }
  
  # plot
  
  CellVol_toplot <- data.frame(CellVol)
  CellVol_toplot$TimePoint <- seq(1,290,1)
  # replace NaN with 0
  CellVol_toplot[is.na(CellVol_toplot)] <- 0
  
  
  pdf(paste0(filepaths.out, "/CellVolume_TimeSeries.pdf"),width=15,height=4,paper='special') 
  print(ggplot(data=melt(CellVol_toplot, id="TimePoint"),
               aes(x=TimePoint, y=value, colour=variable)) +
          geom_line() +
          ylim(0, 35000))
  dev.off()
  
  
  ## Mean growth rate (based on LOWESS and linear fit)
  ############################################################
  
  SpheroidGrowthRate <- list()
  for (pos in 1:38){
    SpheroidGrowthRate[[pos]] <- data.frame(Nr_LOWESS=rep(NA, length = n.rows[[pos]]),
                                            InstantGrowthRate=rep(NA, length = n.rows[[pos]]),
                                            stringsAsFactors=FALSE)
    try(SpheroidGrowthRate[[pos]]$Nr_LOWESS <- lowess(CellVol[1:n.rows[[pos]],pos], f=2/3, iter=3L, delta=0.01 * diff(range(CellVol[1:n.rows[[pos]],pos])))$y, silent = T)
    SpheroidGrowthRate[[pos]]$InstantGrowthRate[-1] <- diff(SpheroidGrowthRate[[pos]]$Nr_LOWESS)
  }
  
  MeanSpheroidGrowthRate <- matrix(nrow=38, ncol=4)
  colnames(MeanSpheroidGrowthRate) <- c("Mean", "STD", "LinearFit_slope", "LinearFit_intercept")
  for (pos in 1:38){
    MeanSpheroidGrowthRate[pos,1] <- mean(SpheroidGrowthRate[[pos]]$InstantGrowthRate[-1])
    MeanSpheroidGrowthRate[pos,2] <- sd(SpheroidGrowthRate[[pos]]$InstantGrowthRate[-1])
    if (sum(!is.na(CellVol[1:n.rows[[pos]],pos]))>0){
      linfit <- lm(CellVol[1:n.rows[[pos]],pos] ~ seq(1,n.rows[[pos]],1))
      MeanSpheroidGrowthRate[pos,3] <- linfit[[1]][2]
      MeanSpheroidGrowthRate[pos,4] <- linfit[[1]][1]
    }
  }
  
  # plot
  
  pdf(paste0(filepaths.out, "/MeanSpheroidGrowthRate.pdf"),width=5,height=4,paper='special') 
  ggplot(data=as.data.frame(MeanSpheroidGrowthRate),
         aes(x=LinearFit_slope)) +
    geom_histogram(colour="black", fill="grey", binwidth = 5) +
    xlim(-3,50)
  dev.off()
  

  ############################################################
  ########## Analysis of nuclei migration           ##########
  ############################################################
  
  NucleiMigrationIn <- list()
  for (pos in 1:38){
    NucleiMigrationIn[[pos]] <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_NucleiMigration/Pos%s.csv", topfolderlocation, PlateIDs[i], pos),
                                           sep = ",", nrows = 290, skip = 1)
    
  }
  
  NucleiMigration <- list()
  for (pos in 1:38){
    NucleiMigration[[pos]] = cbind(NucleiMigrationIn[[pos]][,c(2,3,4)])
    colnames(NucleiMigration[[pos]]) <- c("Time Point", "MeanVelocity", "STDVelocity")
  }
  
  # not all positions have the same number of time points -> determine n.rows for each
  n.rows <- list()
  for (pos in 1:38){
    n.rows[[pos]] <- dim(NucleiMigration[[pos]])[1]
  }
  
  
  ## Mean and STD of track velocity
  ############################################################
  
  FeatureName = c("MeanVelocity", "STDVelocity")
  
  for (fi in 1:2){
    
    FeatureIdx = fi
    
    # user-defined function to label data frame columns by position
    posfunc <- function(x){
      y <- sprintf("Pos_%s", x)
      return(y)
    }
    colnames(SpheroidFeat) <- sapply(X = seq(1,38,1), FUN = posfunc)
    SpheroidFeat <- matrix(NA, ncol=38, nrow=290)  # initialise
    
    for (pos in 1:18){
      SpheroidFeat[1:n.rows[[pos]],pos] <- NucleiMigration[[pos]][,FeatureName[FeatureIdx]]
    }
    for (pos in 21:38){
      SpheroidFeat[1:n.rows[[pos]],pos] <- NucleiMigration[[pos]][,FeatureName[FeatureIdx]]
    }
    
    # plot
    
    SpheroidFeat_toplot <- data.frame(SpheroidFeat)
    SpheroidFeat_toplot$TimePoint <- seq(1,290,1)
    
    pdf(paste0(filepaths.out, "/", FeatureName[FeatureIdx], "_TimeSeries.pdf"),width=15,height=12,paper='special') 
    print(ggplot(data=melt(SpheroidFeat_toplot, id="TimePoint"),
                 aes(x=TimePoint, y=value, colour=variable)) +
            geom_line())
    dev.off()
    
    # violin plot across all time points
    pdf(paste0(filepaths.out, "/", FeatureName[FeatureIdx], "_VlnPlot.pdf"),width=15,height=4,paper='special') 
    print(ggplot(data=melt(SpheroidFeat_toplot, id="TimePoint"),
                 aes(x=variable, y=value, fill=variable)) +
            geom_violin() +
            guides(fill=FALSE) +
            ylim(0,3))
    dev.off()
  }
  

  ############################################################
  ########## Analysis of nuclei distance to center  ##########
  ############################################################
  
  NucleiDistCentIn <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_NucleiDistanceToCenter/MitosisPositioning.csv", topfolderlocation, PlateIDs[i]),
                                 sep = ",", skip=1,
                                 stringsAsFactors = F)[,c(4,3,2,8,11)] 
  colnames(NucleiDistCentIn) <-  c("PosNr", "PosName", "CellPhase","MedianDistance", "STDDistance")
  
  # plot all positions separately
  pdf(paste0(filepaths.out, "/NucleiDistanceToCenter_byCellPhase_byPosition.pdf"),width=15,height=4,paper='special') 
  print(ggplot(data=NucleiDistCentIn,
               aes(x=PosNr, y=MedianDistance, colour=CellPhase)) +
          geom_point(size=2, shape=16) +
          geom_errorbar(aes(ymin=MedianDistance-STDDistance, ymax=MedianDistance+STDDistance), width=.2)
  )
  dev.off()
  
  
  # Function to output summary statistics (median and +/- sd)
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  # violin plot for all positions combined
  pdf(paste0(filepaths.out, "/NucleiDistanceToCenter_byCellPhase_VlnPlot.pdf"),width=5,height=4,paper='special') 
  print(ggplot(data=NucleiDistCentIn,
               aes(x=CellPhase, y=MedianDistance, fill=CellPhase)) +
          geom_violin() +
          stat_summary(fun.data=data_summary) +  # using custom summary function
          guides(fill=FALSE)
  )
  dev.off()
  

  ############################################################
  ########## Analysis of phase transitions          ##########
  ############################################################
  
  
  ## Duration in different phases & nr of transitions 
  ############################################################
  
  PhaseTransitions <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_ClassTransition/ClassDuration.csv", topfolderlocation, PlateIDs[i]),
                                 sep = ",", skip=1)[,c(3,2,4,7,6,5)]
  
  colnames(PhaseTransitions) <-  c("PosName", "CellPhase", "TotalTime",  "NrTransitions", "MedianTimeToTransition", "MaxTimeToTransition") 
  PhaseTransitions$PosNr <- sub("Pos","",PhaseTransitions$PosName)
  
  
  # plot all positions separately
  pdf(paste0(filepaths.out, "/MaxTimeToTransition.pdf"),width=15,height=4,paper='special') 
  ggplot(data=PhaseTransitions,
         aes(x=PosNr, y=MaxTimeToTransition, colour=CellPhase)) +
    geom_point(size=2, shape=16) +
    ylim(0,80)
  dev.off()
  
  
  
  # violin plot for all positions combined
  pdf(paste0(filepaths.out, "/MaxTimeToTransition_VlnPlot.pdf"),width=7,height=4,paper='special') 
  ggplot(data=PhaseTransitions,
         aes(x=CellPhase, y=NrTransitions, fill=CellPhase)) +
    geom_violin() +
    ylim(0,50) +
    stat_summary(fun.data=data_summary) +  # using custom summary function
    theme(axis.text.x = element_blank())
  dev.off()
  
  
  ## T/F transitions 
  ############################################################
  
  PhaseTransitionsTF <- read.table(sprintf("%s/TL%s_S_hSPIM_Output/TimeLapsFiles/Analysis_ClassTransition/Transitions_Quality.csv", topfolderlocation, PlateIDs[i]),
                                   sep = ",", skip=1)[,c(2,3,4)]
  
  colnames(PhaseTransitionsTF) <-  c("PosName", "TF", "count") 
  PhaseTransitionsTF$PosNr <- as.numeric(sub("Pos","",PhaseTransitionsTF$PosName))
  
  pdf(paste0(filepaths.out, "/TF_Transitions.pdf"),width=10,height=6,paper='special') 
  ggplot(data = setorder(PhaseTransitionsTF, "PosNr"),
         aes(x = PosNr, y = count, fill = TF)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous()
  dev.off()
  
  
  ####################################################################
  ##### Collect values in large matrix for downstream analysis #######
  ####################################################################
  
  SummaryMatrix <- matrix(nrow=38, ncol=24)
  
  for (pos in 1:38){
    
    SummaryMatrix[pos,1] <- MeanGrowthRate[pos,1]   # [pos,3] is LOWESS, [pos,1] is average of instantaneous values
    
    SummaryMatrix[pos,2] <- sum(NucleiProperties[[pos]][,4])/sum(NucleiProperties[[pos]][,2])
    SummaryMatrix[pos,3] <- sum(NucleiProperties[[pos]][,5])/sum(NucleiProperties[[pos]][,2])
    SummaryMatrix[pos,4] <- sum(NucleiProperties[[pos]][,6])/sum(NucleiProperties[[pos]][,2])
    
    # Average size of nuclei across all phases
    s <- matrix(nrow=290,ncol=5)
    for (t in 1:290){
      s[t,1] <- NucleiProperties[[pos]][t,3]*NucleiProperties[[pos]][t,7]
      s[t,2] <- NucleiProperties[[pos]][t,4]*NucleiProperties[[pos]][t,8]
      s[t,3] <- NucleiProperties[[pos]][t,5]*NucleiProperties[[pos]][t,9]
      s[t,4] <- NucleiProperties[[pos]][t,6]*NucleiProperties[[pos]][t,10]
      s[t,5] <- sum(s[t,c(1,2,3,4)])/NucleiProperties[[pos]][t,2]
    }      
    SummaryMatrix[pos,5] <- mean(s[,5], na.rm=T)
    
    SummaryMatrix[pos,6] <- mean(NucleiProperties[[pos]][,8], na.rm=T)
    SummaryMatrix[pos,7] <- mean(NucleiProperties[[pos]][,9], na.rm=T)
    SummaryMatrix[pos,8] <- mean(NucleiProperties[[pos]][,10], na.rm=T)
    SummaryMatrix[pos,9] <- mean(NucleiProperties[[pos]][,7], na.rm=T)
    
    SummaryMatrix[pos,10] <- mean(SpheroidProperties[[pos]][,2], na.rm=T)
    
    SummaryMatrix[pos,11] <- mean(CellVol[,pos], na.rm=T)
    
    SummaryMatrix[pos,12] <- mean(MeanSpheroidGrowthRate[pos,3], na.rm=T)  # [pos,3] is LOWESS, [pos,1] is average of instantaneous values
    
    SummaryMatrix[pos,13] <- mean(SpheroidProperties[[pos]][,3], na.rm=T)
    SummaryMatrix[pos,14] <- mean(SpheroidProperties[[pos]][,4], na.rm=T)
    SummaryMatrix[pos,15] <- mean(SpheroidProperties[[pos]][,5], na.rm=T)
    
    SummaryMatrix[pos,16] <- mean(NucleiMigration[[pos]][,2], na.rm=T)
    
    PhaseTransitions <- as.data.table(PhaseTransitions)
    PhaseTransitionsTF <- as.data.table(PhaseTransitionsTF)
    
    SummaryMatrix[pos,17] <- mean(data.matrix(PhaseTransitions[PosNr == pos & CellPhase == "Interphase->Interphase" & MedianTimeToTransition>0 , 5]), na.rm=T)
    SummaryMatrix[pos,18] <- mean(data.matrix(PhaseTransitions[PosNr == pos & CellPhase == "Prophase->Prophase" & MedianTimeToTransition>0 , 5]), na.rm=T)
    SummaryMatrix[pos,19] <- mean(data.matrix(PhaseTransitions[PosNr == pos & CellPhase == "Metphase->Metaphase" & MedianTimeToTransition>0 , 5]), na.rm=T)
    SummaryMatrix[pos,20] <- mean(data.matrix(PhaseTransitions[PosNr == pos & CellPhase == "Anaphase->Anaphase" & MedianTimeToTransition>0 , 5]), na.rm=T)
    
    SummaryMatrix[pos,21] <- sum(data.matrix(PhaseTransitions[PosNr == pos , 4]))
    
    if (length(data.matrix(PhaseTransitionsTF[PosNr == pos]))>0) {    # position exists in table
      if (length(data.matrix(PhaseTransitionsTF[PosNr == pos & TF == "true"])) > 0){      # some are scored 'true'
        SummaryMatrix[pos,22] <- data.matrix(PhaseTransitionsTF[PosNr == pos & TF == "true", 3])/data.matrix(sum(PhaseTransitionsTF[PosNr == pos, 3]))
      } else {    # all are scored 'false'
        SummaryMatrix[pos,22] <- 0
      }    
    }
    
    SummaryMatrix[pos,23] <- SummaryMatrix[pos,10]/SummaryMatrix[pos,15]
    
    SummaryMatrix[pos,24] <- min(DeathTimepoints[,pos], na.rm=T)
    
  }
  
  colnames(SummaryMatrix) <- c("GrowthRate", "ProphaseFrac",	"MetaphaseFrac",	"AnaphaseFrac",	"NucleiSizeAll",	"NucleiSizeP",	"NucleiSizeM",	"NucleiSizeA",	"NucleiSizeI",	"SpheroidVolume",	"CellVolume",	"SpheroidGrowthRate",	"Compactness",	"Convexity",	"Sphericity",	"MigrationSpeed",	"TimeToTransitionI",	"TimeToTransitionP",	"TimeToTransitionM",	"TimeToTransitionA",	"NrTransitions",	"FracFalseTransitions", "Size/Sphericity", "DeathTime")
  
  # save summary matrix
  write.csv(SummaryMatrix, file =paste0(filepaths.out, "/SummaryMatrix.csv"))
  
  
  ####################################################################
  ###############   Determine dead cells    ###########################
  ####################################################################
  
  DeathTimes <- matrix(nrow=1,ncol=38)
  for (pos in 1:38){
    DeathTimes[pos] <- min(DeathTimepoints[,pos], na.rm=T)
  }
  
  DyingCells <- is.finite(DeathTimes)
  
  
  ####################################################################
  ##### Exclude outliers (dead cells) and make heat map  #######
  ####################################################################
  
  
  #SummaryMatrix_sel <- SummaryMatrix[!DyingCells,c(1,2,5,10:11,14:16)]
  SummaryMatrix_sel <- SummaryMatrix[!DyingCells,]
  SummaryMatrix_sel <- SummaryMatrix_sel[,-24]
  
  
  ## Normalise each feature based on mean and SD (could later do across all plates)
  ####################################################################
  
  SummaryMatrix_sel <- scale(SummaryMatrix_sel, center=T, scale=T)
  
  # save summary matrix
  write.csv(SummaryMatrix_sel, file =paste0(filepaths.out, "/SummaryMatrix_exclDyingCells.csv"))
  
  
  # Heatmaps
  
  hm.col <- colorRamps::matlab.like2(101)
  
  # (i) no clustering
  pdf(paste0(filepaths.out,"/24_Factors_Heatmap_noClustering.pdf"),width=20,height=10,paper='special') 
  heatmap.2(SummaryMatrix_sel,
            scale = "none",
            dendrogram = "none",
            na.rm = T,
            Rowv = F,
            Colv = F,
            margins = c(10, 5),
            cexRow = 0.2 + 1/log10(dim(SummaryMatrix_sel)[1]),
            cexCol = 0.2 + 1/log10(dim(SummaryMatrix_sel)[2]),
            col = hm.col
  )
  dev.off()
  
  # (ii) column and row clustering
  pdf(paste0(filepaths.out,"/24_Factors_Heatmap_RowColClustering.pdf"),width=20,height=10,paper='special') 
  heatmap.2(SummaryMatrix_sel,
            scale = "none",
            na.rm = T,
            Rowv = T,
            Colv = T,
            margins = c(10, 5),
            cexRow = 0.2 + 1/log10(dim(SummaryMatrix_sel)[1]),
            cexCol = 0.2 + 1/log10(dim(SummaryMatrix_sel)[2]),
            col = hm.col
  )
  dev.off()
  
  
  
  
  
  
  ####################################################################
  ##### Identify outliers  #######
  ####################################################################
  
  Outliers <- list() # row indices of outliers in SummaryMatrix columns
  
  # outliers = lie outside 1.5*IQR 
  
  for (f in 1:dim(SummaryMatrix)[2]){
    
    M <- median(SummaryMatrix[,f], na.rm = T)
    UpLim <- 3*(quantile(SummaryMatrix[,f], probs = 0.9, na.rm = T) - M) + M
    LoLim <- 3*(quantile(SummaryMatrix[,f], probs = 0.1, na.rm = T) - M) + M
    
    Outliers[[f]] <- c(which(SummaryMatrix[,f] > UpLim), 
                       which(SummaryMatrix[,f] < LoLim)
    )
    
  }
  
  # save outliers
  save(Outliers, file =paste0(filepaths.out, "/Outliers.Rdata"))
  
  
  
  
  
  
  
  ####################################################################
  ##### Check for consistency: correlations  #######
  ####################################################################
  
  ## between the two replicates per target gene
  ####################################################################
  
  # Initialise variables
  linearMod_R2 <- list()
  
  # Replicate IDs
  
  siRNApos_genes_thisPlate_loc <- siRNApos[,siRNAposColIDs[i]]
  siRNApos_genes_thisPlate <- factor(unique(siRNApos_genes_thisPlate_loc[!(siRNApos_genes_thisPlate_loc %in% c("Beads", "empty"))]))
  ReplicateSIRNAs_Idx <-  data.frame(GeneName=siRNApos_genes_thisPlate)
  ReplicateSIRNAs_Idx$ID1 <- match(ReplicateSIRNAs_Idx$GeneName, siRNApos_genes_thisPlate_loc) - 1
  ReplicateSIRNAs_Idx$ID2 <- 41 - match(ReplicateSIRNAs_Idx$GeneName, rev(siRNApos_genes_thisPlate_loc)) - 1
  
  png(filename = paste0(filepaths.out, "/ReplicateCorrelations.png"),
      width = 25, height = 20, units = "cm", res = 300, pointsize = 12)
  par(mfrow=c(4,6))
  ReplicateSIRNAs <- list()
  ReplicateSIRNAs_Correlation <- list()
  for (f in 1:22){   # 22 being the number of factors
    ReplicateSIRNAs[[f]] <- matrix(ncol=2, nrow=dim(ReplicateSIRNAs_Idx)[1])
    for (j in 1:dim(ReplicateSIRNAs_Idx)[1]){
      ReplicateSIRNAs[[f]][j,1] <- SummaryMatrix[ReplicateSIRNAs_Idx$ID1[j],f]
      ReplicateSIRNAs[[f]][j,2] <- SummaryMatrix[ReplicateSIRNAs_Idx$ID2[j],f]
    }
    
    ReplicateSIRNAs_Correlation[[f]] <- cor(ReplicateSIRNAs[[f]][,1],ReplicateSIRNAs[[f]][,2])
    
    ReplicateSIRNAs_Regression.df <- data.frame(SampleA = ReplicateSIRNAs[[f]][,1],
                                                SampleB = ReplicateSIRNAs[[f]][,2])
    
    linearMod <- lm(SampleA ~ SampleB, data=ReplicateSIRNAs_Regression.df)
    # print(linearMod)
    # summary(linearMod) 
    linearMod_R2[[f]] <- summary(lm(SampleA ~ SampleB, ReplicateSIRNAs_Regression.df))$r.squared
    
    plot(ReplicateSIRNAs[[f]], 
         main = paste0(colnames(SummaryMatrix)[f], "\nR = ", round(linearMod_R2[[f]], digits = 2))
    )
    
  }
  dev.off()
  
  
  ## plot without outliers
  ####################################################################
  
  linearMod_noOutliers_R2 <- list()
  
  png(filename = paste0(filepaths.out, "/ReplicateCorrelations_noOutliers.png"),
      width = 25, height = 20, units = "cm", res = 300, pointsize = 12)
  
  par(mfrow=c(4,6))
  for (f in 1:22){   # 22 being the number of factors
    OutlierRows <- as.logical((ReplicateSIRNAs_Idx$ID1 %in% Outliers[[f]])+((ReplicateSIRNAs_Idx$ID2 %in% Outliers[[f]])>0))
    
    ReplicateSIRNAs_noOutliers_Regression.df <- data.frame(SampleA = ReplicateSIRNAs[[f]][!OutlierRows,1],
                                                           SampleB = ReplicateSIRNAs[[f]][!OutlierRows,2])
    
    linearMod_noOutliers <- lm(SampleA ~ SampleB, data=ReplicateSIRNAs_noOutliers_Regression.df)
    # print(linearMod)
    # summary(linearMod) 
    linearMod_noOutliers_R2[[f]] <- summary(lm(SampleA ~ SampleB, ReplicateSIRNAs_noOutliers_Regression.df))$r.squared
    
    plot(ReplicateSIRNAs[[f]][!OutlierRows,], 
         main = paste0(colnames(SummaryMatrix)[f], "\nR = ", round(linearMod_noOutliers_R2[[f]], digits = 2))
    )
    # xlim = c(0,20), ylim=c(0,20))
    
  }
  
  dev.off()
  
  
  
  
  
}


####################################################################
####################################################################
####################################################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
