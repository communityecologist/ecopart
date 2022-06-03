#' Partitioning pairwise dissimilarity
#'
#' \code{ecopart.pair} patitions the temporal changes in pairwise dissimilarity (Jaccard, Sorensen, Ruzicka, and Bray-Curtis indices) into dynamic components based on methods proposed by Tatsumi et al. (2021, 2022)
#'
#' @param d1 A matrix or dataframe at time 1. Rows are a pair of sites (sites 1 and 2), columns are species, and elements are presence-absence (01) or abundances of species.
#' @param d2 A matrix or dataframe at time 2. Note that \code{d1} and \code{d2} must have exactly the same sites and species in the same order.
#' @param index Type of dissimilarity measure. Options are \code{"jaccard"}, \code{"sorensen"}, \code{"ruzicka"}, and \code{"bray-curtis"}.
#' \describe{
#'   \item{\code{"jaccard"} and \code{"sorensen"}}{These indices are based on presence-absence data (Jaccard 1912, Sorensen 1948). When \code{d1} and \code{d2} are abundance data, the elements are automatically converted to presence-absence data by replacing non-zero values with 1.}
#'   \item{\code{"ruzicka"} and \code{"bray-curtis"}}{These indices are based on abundance data (Bray & Curtis 1957, Ruzicka 1958). When \code{d1} and \code{d2} are presence-absence data, Ruzicka and Bray-Curtis indices are equivalent to the Jaccard and Sorensen indices, respectively.}
#' }
#' @param components Types of components into which the total change in beta diversity is partitioned. Options are \code{"two"}, \code{"four"}, \code{"six"}, and \code{"sp"}.
#' \describe{
#'   \item{\code{"two"}}{Calculates extinction and colonization components (when \code{index = "jaccard"} or \code{index = "sorensen"}) or subtractive and additive components (when \code{"ruzicka"} or \code{"bray-curtis"}).}
#'   \item{\code{"four"}}{Calculates extinction- and colonization-induced homogenization and differentiation (when \code{index = "jaccard"} or \code{index = "sorensen"}) or subtractive and additive homogenization and differentiation (when \code{"ruzicka"} or \code{"bray-curtis"}).}
#'   \item{\code{"six"}}{Calculates six components that represent (1) homogenization by decrease in U, (2) differentiation by decrease in C, (3) differentiation by change from C to U, (4) differentiation by increase in U, (5) homogenization by increase in C, and (6) homogenization by change from U to C, where U is the existance or abundance of species unique to either site and C is the existance or abundance of species common to both sites.}
#'   \item{\code{"sp"}}{Same with \code{"four"} but the components are further partitioned down to the species-level.}
#' }
#' @return The \code{ecopart.pair} function returns a vector or matrix object containing the partitioned components of beta diversity.
#' * When \code{components = "two"} was specified, the function returns a vector object with two elements: extinction and colonization components (when \code{index = "jaccard"} or \code{index = "sorensen"}) or subtractive and additive components (when \code{index = "ruzicka"} or \code{index = "bray-curtis"}). The extinction and colonization components represent temporal changes in beta diversity that result from local species extinctions and colonizations. The subtractive and additive components represent temporal changes in beta diversity that result from local losses and gains in species abundances.
#' * When \code{components = "four"} was specified, the function returns a vector object with four elements: extinction- and colonization-induced homogenization and differentiation (when \code{index = "jaccard"} or \code{index = "sorensen"}) or subtractive and additive homogenization and differentiation (when \code{index = "ruzicka"} or \code{index = "bray-curtis"}). Homogenization and differentiation indicate decreases and increases in beta diversity, respectively.
#' * When \code{components = "six"} was specified, the function returns a vector object with six elements: (1) homogenization by decrease in U, (2) differentiation by decrease in C, (3) differentiation by change from C to U, (4) differentiation by increase in U, (5) homogenization by increase in C, and (6) homogenization by change from U to C.
#' * When \code{components = "sp"} was specified, the function returns a matrix object. The rows represent the four components that are equivalent to when \code{components = "four"} was specified. The columns represent species.
#' @author Shinichi Tatsumi
#' @references
#' * Bray JR, Curtis JT (1957) An ordination of the upland forest communities of southern Wisconsin. \emph{Ecological Monographs} 27(4): 325-349.
#' * Jaccard P (1912) The distribution of the flora in the alphine zone. \emph{New Phytologist} 11(2): 37-50.
#' * Ruzicka M (1958) Anwendung mathematisch-statisticher Methoden in der Geobotanik (synthetische Bearbeitung von Aufnahmen). \emph{Biologia, Bratislava} 13: 647-661.
#' * Sorensen T (1948) A method of establishing groups of equal amplitude in plant sociology based on similarity of species and its application to analyses of the vegetation on Danish commons. \emph{Kongelige Danske Videnskabernes Selskab} 5(4): 1-34.
#' * Tatsumi S, Iritani R, Cadotte MW (2021) Temporal changes in spatial variation: partitioning the extinction and colonisation components of beta diversity. \emph{Ecology Letters} 24(5): 1063-1072.
#' * Tatsumi S, Iritani R, Cadotte MW (2022) Partitioning the temporal changes in abundance-based beta diversity into loss and gain components. \emph{Methods in Ecology and Evolution}: in press
#' @export

ecopart.pair  <- function(d1, d2, index="sorensen", components="four"){

  components <-	match.arg(components, c("two", "four", "six", "sp"))
  index      <-	match.arg(index, c("jaccard", "sorensen", "ruzicka", "bray-curtis"))
  N          <-	nrow(d1)
  S          <-	ncol(d1)

  if((index=="jaccard")|(index=="sorensen")){
    d1[d1!=0] <- 1
    d2[d2!=0] <- 1
  }

  Mat   <- rbind(d1, d2)
  Ord   <- apply(Mat, 2, FUN="order")
  Srt   <- apply(Mat, 2, FUN="sort")

  LossL <- ifelse(colSums(Ord[1:2,])==7 | colSums(Ord[1:2,])!=3 & Ord[3,]>Ord[4,], Srt[4,]-Srt[3,], 0)
  LossE <- ifelse(colSums(Ord[1:2,])==7,                                           Srt[3,]-Srt[2,], 0)
  LossS <- ifelse(colSums(Ord[1:2,])==7 | colSums(Ord[1:2,])!=3 & Ord[1,]>Ord[2,], Srt[2,]-Srt[1,], 0)
  GainL <- ifelse(colSums(Ord[1:2,])==3 | colSums(Ord[1:2,])!=7 & Ord[3,]<Ord[4,], Srt[4,]-Srt[3,], 0)
  GainE <- ifelse(colSums(Ord[1:2,])==3,                                           Srt[3,]-Srt[2,], 0)
  GainS <- ifelse(colSums(Ord[1:2,])==3 | colSums(Ord[1:2,])!=7 & Ord[1,]<Ord[2,], Srt[2,]-Srt[1,], 0)
  HddnD <- ifelse(colSums(Ord[1:2,])==5,                                           Srt[3,]-Srt[2,], 0)

  C1  <- apply(d1, 2, min) |> sum()
  C2  <- apply(d2, 2, min) |> sum()
  U1  <- as.matrix(d1) |> diff() |> abs() |> sum()
  U2  <- as.matrix(d2) |> diff() |> abs() |> sum()
  X   <- ifelse(index=="jaccard"|index=="ruzicka", C1*U1/(C1+U1)/(C2+U2), 2*C1*U1/(2*C1+U1)/(2*C2+U2))

  DBeta     <-  matrix(nrow=6, ncol=S)
  DBeta[1,]	<- -X/U1          * (LossL + HddnD)
  DBeta[2,]	<-  X/C1          *  LossE
  DBeta[3,]	<- (X/C1 + X/U1)  *  LossS
  DBeta[4,]	<-  X/U1          * (GainL + HddnD)
  DBeta[5,]	<- -X/C1          *  GainE
  DBeta[6,]	<- -(X/C1 + X/U1) *  GainS

  if(components=="two"){

    Res          <-  c(sum(DBeta[1:3,]), sum(DBeta[4:6,]))
    if(index=="jaccard"|index=="sorensen"){
      names(Res) <-  c("Extinction component", "Colonization component")
    } else if(index=="ruzicka"|index=="bray-curtis"){
      names(Res) <-  c("Subtractive component", "Additive component")
    }

  } else if(components=="four"){

    Res          <-  c(sum(DBeta[1,]), sum(DBeta[2:3,]), sum(DBeta[5:6,]), sum(DBeta[4,]))
    if(index=="jaccard"|index=="sorensen"){
      names(Res) <-  c("Extinction homogenization", "Extinction differentiation",
                       "Colonization homogenization", "Colonization differentiation")
    } else if(index=="ruzicka"|index=="bray-curtis"){
      names(Res) <-  c("Subtractive homogenization", "Subtractive differentiation",
                       "Additive homogenization", "Additive differentiation")
    }

  } else if(components=="six"){

    Res          <-  rowSums(DBeta)
    names(Res)   <-  c("DeltaBeta1", "DeltaBeta2", "DeltaBeta3", "DeltaBeta4", "DeltaBeta5", "DeltaBeta6")

  } else if(components=="sp"){

    Res             <-  sapply(1:ncol(DBeta), function(s) c(DBeta[1,s], sum(DBeta[2:3,s]), sum(DBeta[5:6,s]), DBeta[4,s]))
    if(index=="jaccard"|index=="sorensen"){
      rownames(Res) <-  c("Extinction homogenization", "Extinction differentiation",
                          "Colonization homogenization", "Colonization differentiation")
    } else if(index=="ruzicka"|index=="bray-curtis"){
      rownames(Res) <-  c("Subtractive homogenization", "Subtractive differentiation",
                          "Additive homogenization", "Additive differentiation")
    }
    colnames(Res)   <-  colnames(d1)

  }

  return(Res)
}
