#####################################################
##   R functions for extracting the Swedish fish   ##
##   community data from RivFishTIME               ##
##   by Shinichi Tatsumi                           ##
##                             December 29, 2021   ##
#####################################################

# This code is used in:
# Tatsumi S, Iritani R, Cadotte MW (2021) Partitioning the temporal changes in abundance-based beta diversity into loss and gain components. bioRxiv, 2021.05.28.446095.

SiteInfo    <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv")
Dat         <- read.csv("1873_2_RivFishTIME_SurveyTable.csv")

SiteInfo           <- SiteInfo[SiteInfo$SourceID==42,]
SiteInfo$Region    <- factor(SiteInfo$Region, levels=sort(unique(SiteInfo$Region)))
SiteInfo$Province  <- factor(SiteInfo$Province, levels=sort(unique(SiteInfo$Province)))
SiteInfo$Waterbody <- factor(SiteInfo$Waterbody, levels=sort(unique(SiteInfo$Waterbody)))
Dat         <- Dat[!is.na(match(Dat$TimeSeriesID, SiteInfo$TimeSeriesID)),]
Dat         <- cbind(Dat, SiteInfo[match(Dat$TimeSeriesID, SiteInfo$TimeSeriesID),])
Dat         <- Dat[Dat$Quarter==3 | Dat$Quarter==4,]
Dat$Quarter <- factor(Dat$Quarter, levels=sort(unique(Dat$Quarter)))
Dat$Species <- factor(Dat$Species, levels=sort(unique(Dat$Species)))
Site        <- as.character(sort(unique(Dat$SiteID)))
Year        <- sort(unique(Dat$Year))
SiteYear    <- unique(paste(Dat$SiteID, Dat$Year, sep="_"))
Survey      <- data.frame(matrix(unlist(strsplit(SiteYear, "_")), byrow=T, ncol=2))
Survey[,1]  <- as.character(Survey[,1])
Survey[,2]  <- as.numeric(as.character(Survey[,2]))
Mat         <- matrix(0, nrow=length(Site), ncol=length(Year))
rownames(Mat) <- Site
colnames(Mat) <- Year
for(i in 1:nrow(Survey)){
  Temp <- sum(as.numeric(unique(Dat[which(Dat$SiteID==Survey[i,1] & Dat$Year==Survey[i,2]), "Quarter"])))
  Mat[match(Survey[i,1], Site), match(Survey[i,2], Year)] <- Temp # Q3=1, Q4=2, Q3&Q4=3
}
Y1 <- 1990; Y2 <- 2018
Mat6      <- Mat[,Year%in%c(Y1, Y2)]
Region    <- (SiteInfo[match(Site, SiteInfo$SiteID), "Region"])[sapply(1:nrow(Mat6), function(x) all(Mat6[x,]>0))]
Province  <- (SiteInfo[match(Site, SiteInfo$SiteID), "Province"])[sapply(1:nrow(Mat6), function(x) all(Mat6[x,]>0))]
Waterbody <- (SiteInfo[match(Site, SiteInfo$SiteID), "Waterbody"])[sapply(1:nrow(Mat6), function(x) all(Mat6[x,]>0))]
Mat6      <- Mat6[sapply(1:nrow(Mat6), function(x) all(Mat6[x,]>0)),]

RPW        <- paste(Region, Province, Waterbody, sep="_")
WBselected <- names(xtabs(~RPW)[xtabs(~RPW)>=2])
Mat6       <- Mat6[!is.na(match(RPW, WBselected)),]
RPW        <- RPW[!is.na(match(RPW, WBselected))]
D1 <- matrix(0, nrow=nrow(Mat6), ncol=length(unique(Dat$Species)))
D2 <- matrix(0, nrow=nrow(Mat6), ncol=length(unique(Dat$Species)))
rownames(D1) <- rownames(Mat6); colnames(D1) <- sort(unique(Dat$Species))
rownames(D2) <- rownames(Mat6); colnames(D2) <- sort(unique(Dat$Species))

for(i in 1:nrow(Mat6)){
  Q    <- ifelse(Mat6[i,1]==1, 3, 4)
  Temp <- Dat[Dat$Year==Y1 & Dat$SiteID==rownames(Mat6)[i] & Dat$Quarter==Q, c("Species", "Abundance")]
  for(j in 1:nrow(Temp)){
    D1[i, which(colnames(D1)==Temp$Species[j])] <- Temp$Abundance[j]
  }
  Q    <- ifelse(Mat6[i,2]==1, 3, 4)
  Temp <- Dat[Dat$Year==Y2 & Dat$SiteID==rownames(Mat6)[i] & Dat$Quarter==Q, c("Species", "Abundance")]
  for(j in 1:nrow(Temp)){
    D2[i, which(colnames(D2)==Temp$Species[j])] <- Temp$Abundance[j]
  }
}

# 01 Harrison-Baselga
Dpa1 <- D1
Dpa1[Dpa1>0] <- 1
Dpa2 <- D2
Dpa2[Dpa2>0] <- 1

aaa <-  sapply(1:length(unique(RPW)), function(w) hb(Dpa1[RPW==unique(RPW)[w],]))
bbb <-  sapply(1:length(unique(RPW)), function(w) hb(Dpa2[RPW==unique(RPW)[w],]))
XXX <-  bbb-aaa
HBs <- sapply(1:length(unique(RPW)), function(w) ecopart.multi.abund(Dpa1[RPW==unique(RPW)[w],], Dpa2[RPW==unique(RPW)[w],]))

# Abundance Harrison-Baselga
aaa <-  sapply(1:length(unique(RPW)), function(w) hb(D1[RPW==unique(RPW)[w],]))
bbb <-  sapply(1:length(unique(RPW)), function(w) hb(D2[RPW==unique(RPW)[w],]))
XXX <-  bbb-aaa
HBs <- sapply(1:length(unique(RPW)), function(w) ecopart.multi.abund(D1[RPW==unique(RPW)[w],], D2[RPW==unique(RPW)[w],]))
