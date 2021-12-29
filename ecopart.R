#########################################################
##   R functions ecopart.pair() and ecopart.multi()    ##
##   by Shinichi Tatsumi                               ##
##                                December 29, 2021    ##
#########################################################

# CC-BY-4.0
# Tatsumi S, Iritani R, Cadotte MW (2021) Temporal changes in spatial variation: partitioning the extinction and colonisation components of beta diversity. Ecology Letters 24(5): 1063???1072.
# Tatsumi S, Iritani R, Cadotte MW (2021) Partitioning the temporal changes in abundance-based beta diversity into loss and gain components. bioRxiv, 2021.05.28.446095.

ecopart.pair	<-	function(d1, d2, index="ruzicka", components="two"){
  
  components <-	match.arg(components, c("two", "four", "six", "sp"))
  index		   <-	match.arg(index, c("ruzicka", "bray-curtis", "harrison-baselga"))
  N				   <-	nrow(d1)
  S			     <-	ncol(d1)
  p          <- ifelse(index=="ruzicka", 1, 2)
  
  if(components=="two"){
    Res        <- lapply(1:2, function(x) matrix(nrow=N, ncol=N))
    names(Res) <- c("Changes in beta by abundance loss (ΔBetaLoss)",
                    "Changes in beta by abundance gain (ΔBetaGain)")
  } else if(components=="four"){
    Res        <- lapply(1:4, function(x) matrix(nrow=N, ncol=N))
    names(Res) <- c("Homogenization by abundance loss (ΔBetaLoss-)",
                    "Heterogenization by abundance loss (ΔBetaLoss+)",
                    "Homogenization by abundance gain (ΔBetaGain-)",
                    "Heterogenization by abundance gain (ΔBetaGain+)")
  } else if(components=="six"){
    Res        <- lapply(1:6, function(x) matrix(nrow=N, ncol=N))
    names(Res) <- c("ΔBeta1", "ΔBeta2", "ΔBeta3", "ΔBeta4", "ΔBeta5", "ΔBeta6")
  } else if(components=="sp"){
    Res        <- lapply(1:S, function(x) matrix(nrow=N, ncol=N))
    names(Res) <- colnames(d1)
  }
  
  for(j in 1:(N-1)){
    for(i in (j+1):N){
      
      p1    <- d1[c(i,j),]
      p2    <- d2[c(i,j),]
      Mat	  <- rbind(p1, p2)
      Ord   <- apply(Mat, 2, FUN="order")
      Srt   <- apply(Mat, 2, FUN="sort")
      
      LossL <- ifelse(colSums(Ord[1:2,])==7 | colSums(Ord[1:2,])!=3 & Ord[3,]>Ord[4,], Srt[4,]-Srt[3,], 0)
      LossE <- ifelse(colSums(Ord[1:2,])==7,                                           Srt[3,]-Srt[2,], 0)
      LossS <- ifelse(colSums(Ord[1:2,])==7 | colSums(Ord[1:2,])!=3 & Ord[1,]>Ord[2,], Srt[2,]-Srt[1,], 0)
      GainL <- ifelse(colSums(Ord[1:2,])==3 | colSums(Ord[1:2,])!=7 & Ord[3,]<Ord[4,], Srt[4,]-Srt[3,], 0)
      GainE <- ifelse(colSums(Ord[1:2,])==3,                                           Srt[3,]-Srt[2,], 0)
      GainS <- ifelse(colSums(Ord[1:2,])==3 | colSums(Ord[1:2,])!=7 & Ord[1,]<Ord[2,], Srt[2,]-Srt[1,], 0)
      HddnD <- ifelse(colSums(Ord[1:2,])==5,                                           Srt[3,]-Srt[2,], 0)
      
      if((index=="ruzicka")|(index=="bray-curtis")){
        
        C1  <- sum(apply(p1, 2, min))
        C2  <- sum(apply(p2, 2, min))
        U1  <- sum(abs(diff(p1)))
        U2  <- sum(abs(diff(p2)))
        X   <- p*C1*U1 / (p*C1+U1) / (p*C2+U2)
        
        DBeta			<-　 matrix(nrow=6, ncol=S)
        DBeta[1,]	<-　 -X/U1				 * (LossL + HddnD)
        DBeta[2,]	<-　  X/C1				 *  LossE
        DBeta[3,]	<- 　(X/C1 + X/U1) *  LossS
        DBeta[4,]	<-  　X/U1				 * (GainL + HddnD)
        DBeta[5,]	<- 　-X/C1				 *  GainE
        DBeta[6,]	<-  -(X/C1 + X/U1) *  GainS
      
      } else if(index=="harrison-baselga"){
        
        A  <- sum(p1)
        M  <- sum(apply(p1, 2, max))
        Y  <- 2*M / sum(p2)
        
        DBeta			<-　 matrix(nrow=6, ncol=S)
        DBeta[1,]	<-　 (-Y/M +   Y/A)	* (LossL + HddnD)
        DBeta[2,]	<-　 (-Y/M + 2*Y/A)	*  LossE
        DBeta[3,]	<- 　          Y/A  *  LossS
        DBeta[4,]	<-    (Y/M -   Y/A)	* (GainL + HddnD)
        DBeta[5,]	<- 　 (Y/M - 2*Y/A)	*  GainE
        DBeta[6,]	<-　          -Y/A  *  GainS
      }
      
      if(components=="two"){
        Res[[1]][i,j]  <-  sum(DBeta[1:3,])
        Res[[2]][i,j]  <-  sum(DBeta[4:6,])
      } else if(components=="four"){
        Res[[1]][i,j]  <-  sum(DBeta[1,])
        Res[[2]][i,j]  <-  sum(DBeta[2:3,])
        Res[[3]][i,j]  <-  sum(DBeta[5:6,])
        Res[[4]][i,j]  <-  sum(DBeta[4,])
      } else if(components=="six"){
        for(x in 1:6) Res[[x]][i,j] <-  sum(DBeta[x,])
      } else if(components=="sp"){
        for(x in 1:S) Res[[x]][i,j] <-  sum(DBeta[,x])
      }
      
    } # End loop i
  } # End loop j
  
  for(r in 1:length(Res)) Res[[r]] <- as.dist(Res[[r]])
  return(Res)
} # ecopart.pair.adund()

ecopart.multi <-	function(d1, d2, components="two"){
  
  components <-  match.arg(components, c("two", "four", "sp"))
  N		  	   <-  nrow(d1)
  S			     <-  ncol(d1)
  not.zero   <-  function(x) sapply(1:length(x), function(i) all.equal(as.numeric(x[i]), 0)) != TRUE
  DBeta      <-  matrix(nrow=4, ncol=S)
  
  A  <-  sum(d1)
  M  <-  sum(apply(d1, 2, max))
  Y  <-  N*M / (N-1) / sum(d2)
  
  for(s in 1:S){
    
    Abn1   <-  d1[,s]
    Abn2   <-  d2[,s]
    Max1   <-  max(Abn1)
    Max2   <-  max(Abn2)
    A2_M1  <-  max(Abn2[Abn1==Max1])
    A1_M2  <-  max(Abn1[Abn2==Max2])
    AX     <-  Abn1>A2_M1 & Abn1>A1_M2 & Abn2>A2_M1 & Abn2>A1_M2
    AX     <-  ifelse(any(AX), max(pmin(Abn1[AX], Abn2[AX])), 0)
    LossL  <-  ifelse(Abn1>Max2, Abn1-Max2, 0)
    LossD  <-  ifelse(Max1>A2_M1 & Max2>A1_M2 & Abn1>A2_M1 & Abn1>A1_M2 & Abn1>AX & Abn1>Abn2,
                      pmin(Max1, Max2, Abn1) - max(A2_M1, A1_M2, AX), 0)
    LossS  <-  ifelse(Abn1>Abn2, Abn1-Abn2-LossL-LossD, 0)
    GainL  <-  ifelse(Abn2>Max1, Abn2-Max1, 0)
    GainD  <-  ifelse(Max1>A2_M1 & Max2>A1_M2 & Abn2>A2_M1 & Abn2>A1_M2 & Abn2>AX & Abn2>Abn1,
                      pmin(Max1, Max2, Abn2) - max(A2_M1, A1_M2, AX), 0)
    GainS  <-  ifelse(Abn2>Abn1, Abn2-Abn1-GainL-GainD, 0)
    
    LossL.up   <-  ifelse(not.zero(LossL), Abn1, 0)
    LossL.low  <-  ifelse(not.zero(LossL), Abn1-LossL, 0)
    LossD.up   <-  ifelse(not.zero(LossD), Abn1-LossL, 0)
    LossD.low  <-  ifelse(not.zero(LossD), Abn1-LossL-LossD, 0)
    LossS.up   <-  ifelse(not.zero(LossS), Abn1-LossL-LossD, 0)
    LossS.low  <-  ifelse(not.zero(LossS), Abn1-LossL-LossD-LossS, 0)
    GainL.up   <-  ifelse(not.zero(GainL), Abn2, 0)
    GainL.low  <-  ifelse(not.zero(GainL), Abn2-GainL, 0)
    GainD.up   <-  ifelse(not.zero(GainD), Abn2-GainL, 0)
    GainD.low  <-  ifelse(not.zero(GainD), Abn2-GainL-GainD, 0)
    GainS.up   <-  ifelse(not.zero(GainS), Abn2-GainL-GainD, 0)
    GainS.low  <-  ifelse(not.zero(GainS), Abn2-GainL-GainD-GainS, 0)
    
    Lvl        <-  sort(c(0, Abn1, Abn2))
    Lvl.mean   <-  sapply(1:(2*N), function(i) sum(Lvl[i:(i+1)])/2)
    Lvl.diff   <-  diff(Lvl)
    LossL.n    <-  sapply(1:(2*N), function(i) sum(Lvl.mean[i]<LossL.up & Lvl.mean[i]>LossL.low))
    LossD.n    <-  sapply(1:(2*N), function(i) sum(Lvl.mean[i]<LossD.up & Lvl.mean[i]>LossD.low))
    LossS.n    <-  sapply(1:(2*N), function(i) sum(Lvl.mean[i]<LossS.up & Lvl.mean[i]>LossS.low))
    GainL.n    <-  sapply(1:(2*N), function(i) sum(Lvl.mean[i]<GainL.up & Lvl.mean[i]>GainL.low))
    GainD.n    <-  sapply(1:(2*N), function(i) sum(Lvl.mean[i]<GainD.up & Lvl.mean[i]>GainD.low))
    GainS.n    <-  sapply(1:(2*N), function(i) sum(Lvl.mean[i]<GainS.up & Lvl.mean[i]>GainS.low))
    
    Term.1     <-  ifelse(LossL.n!=0,  Y * (-1/M + LossL.n/A) * Lvl.diff, 0)
    Term.2     <-  ifelse(LossD.n!=0,  Y * (-1/M + LossD.n/A) * Lvl.diff, 0)
    Term.3     <-  ifelse(LossS.n!=0,  Y *         LossS.n/A  * Lvl.diff, 0)
    Term.4     <-  ifelse(GainL.n!=0,  Y * ( 1/M - GainL.n/A) * Lvl.diff, 0)
    Term.5     <-  ifelse(GainD.n!=0,  Y * ( 1/M - GainD.n/A) * Lvl.diff, 0)
    Term.6     <-  ifelse(GainS.n!=0, -Y *         GainS.n/A  * Lvl.diff, 0)
    
    DBeta[1,s] <-  sum(Term.1[Term.1<0], Term.2[Term.2<0])
    DBeta[2,s] <-  sum(Term.1[Term.1>0], Term.2[Term.2>0], Term.3)
    DBeta[3,s] <-  sum(Term.4[Term.4<0], Term.5[Term.5<0], Term.6)
    DBeta[4,s] <-  sum(Term.4[Term.4>0], Term.5[Term.5>0])
  } # End loop s (species)
  
  if(components=="two"){
    Res				    <-	c(sum(DBeta[1:2,]), sum(DBeta[3:4,]))
    names(Res)    <-  c("Changes in beta by abundance loss (ΔBetaLoss)",
                        "Changes in beta by abundance gain (ΔBetaGain)")
  
  } else if(components=="four"){
    Res			      <-	c(sum(DBeta[1,]), sum(DBeta[2,]), sum(DBeta[3,]), sum(DBeta[4,]))
    names(Res)    <-  c("Homogenization by abundance loss (ΔBetaLoss-)",
                        "Heterogenization by abundance loss (ΔBetaLoss+)",
                        "Homogenization by abundance gain (ΔBetaGain-)",
                        "Heterogenization by abundance gain (ΔBetaGain+)")
  
  } else if(components=="sp"){
    Res						<-	DBeta
    rownames(Res)	<-	c("ΔBetaLoss-", "ΔBetaLoss+", "ΔBetaGain-", "ΔBetaGain+")
    colnames(Res)	<-	colnames(d1)
  }
  
  return(Res)
} # ecopart.multi.abund()
