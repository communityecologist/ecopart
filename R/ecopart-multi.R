#' Partitioning multi-site beta diversity
#'
#' \code{ecopart.multi} patitions the temporal changes in multi-site beta diversity (Whittaker or Baselga's beta) into dynamic components based on methods proposed by Tatsumi et al. (2021, 2022)
#'
#' @param d1 Matrix or dataframe at time 1. Rows are sites, columns are species, and elements are presence-absence (01) or abundances of species.
#' @param d2 Matrix or dataframe at time 2. Note that \code{d1} and \code{d2} must have exactly the same sites and species in the same order.
#' @param index Type of beta diversity measure. Options are \code{"whittaker"} and \code{"baselga"}.
#' \describe{
#'   \item{\code{"whittaker"}}{This index is based on presence-absence data (Whittaker 1960). When \code{d1} and \code{d2} are abundance data, the elements are automatically converted to presence-absence data by replacing non-zero values with 1.}
#'   \item{\code{"baselga"}}{This index, also known as the normalized abundance-based Whittaker's beta, is based on abundance data (Baselga 2017). When \code{d1} and \code{d2} are presence-absence data, Baselga's beta is equivalent to Harrison's beta (Harrison et al. 1992), also known as the normalized Whittaker's beta.}
#' }
#' @param components Types of components into which the total change in beta diversity is partitioned. Options are \code{"two"}, \code{"four"}, and \code{"sp"}.
#' \describe{
#'   \item{\code{"two"}}{Calculates extinction and colonization components (when \code{index = "whittaker"}) or subtractive and additive components (when \code{index = "baselga"}).}
#'   \item{\code{"four"}}{Calculates extinction- and colonization-induced homogenization and differentiation (when \code{index = "whittaker"}) or subtractive and additive homogenization and differentiation (when \code{index = "baselga"}).}
#'   \item{\code{"sp"}}{Same with \code{"four"} but the components are further partitioned down to the species-level.}
#' }
#' @return The \code{ecopart.multi} function returns a vector or matrix object containing the partitioned components of beta diversity.
#' * When \code{components = "two"} was specified, the function returns a vector object with two elements: extinction and colonization components (when \code{index = "whittaker"}) or subtractive and additive components (when \code{index = "baselga"}). The extinction and colonization components represent temporal changes in beta diversity that result from local species extinctions and colonizations. The subtractive and additive components represent temporal changes in beta diversity that result from local losses and gains in species abundances.
#' * When \code{components = "four"} was specified, the function returns a vector object with four elements: extinction- and colonization-induced homogenization and differentiation (when \code{index = "whittaker"}) or subtractive and additive homogenization and differentiation (when \code{index = "baselga"}). Homogenization and differentiation indicate decreases and increases in beta diversity, respectively.
#' * When \code{components = "sp"} was specified, the function returns a matrix object. The rows represent the four components that are equivalent to when \code{components = "four"} was specified. The columns represent species.
#' @author Shinichi Tatsumi
#' @references
#' * Baselga A (2017) Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. \emph{Methods in Ecology and Evolution} 8(7): 799-808.
#' * Harrison S, Ross SJ, Lawton JH (1992) Beta diversity on geographic gradients in Britain. \emph{Journal of Animal Ecology} 61(1): 151-158.
#' * Tatsumi S, Iritani R, Cadotte MW (2021) Temporal changes in spatial variation: partitioning the extinction and colonisation components of beta diversity. \emph{Ecology Letters} 24(5): 1063-1072.
#' * Tatsumi S, Iritani R, Cadotte MW (2022) Partitioning the temporal changes in abundance-based beta diversity into loss and gain components. \emph{Methods in Ecology and Evolution}: in press
#' * Whittaker RH (1960) Vegetation of the Siskiyou Mountains, Oregon and California. \emph{Ecological Monographs} 30(3): 279-338.
#' @export

ecopart.multi  <-  function(d1, d2, index="whittaker", components="four"){

  index      <-  match.arg(index, c("whittaker", "baselga"))
  components <-  match.arg(components, c("two", "four", "sp"))
  N          <-  nrow(d1)
  S          <-  ncol(d1)

  if(index=="whittaker"){

    d1[d1!=0] <- 1
    d2[d2!=0] <- 1

    X       <-  colSums(d1)
    Y       <-  colSums(d2)
    Z       <-  sapply(1:S, function(s) sum(d1[,s]==1 & d2[,s]==1))
    Alpha1  <-  rowSums(d1) |> mean()
    Alpha2  <-  rowSums(d2) |> mean()
    Gamma1  <-  sum(colSums(d1) > 0)
    Beta1   <-  Gamma1/Alpha1
    H       <-  Beta1/((Alpha2-Alpha1)/Alpha1 + 1)

    Term.1  <-  (X    *H/N/Alpha1 - H/Gamma1)   * (X>0 & Z==0)
    Term.2  <-  (X-Z) *H/N/Alpha1               * (X>Z & Z>0)
    Term.3  <-  (-Y   *H/N/Alpha1 + H/Gamma1)   * (Y>0 & Z==0)
    Term.4  <-  (Z-Y) *H/N/Alpha1               * (Y>Z & Z>0)

    DBeta       <-  matrix(nrow=4, ncol=S)
    DBeta[1,]   <-  ifelse(Term.1<0, Term.1, 0)
    DBeta[2,]   <-  ifelse(Term.1>0, Term.1, 0) + Term.2
    DBeta[3,]   <-  ifelse(Term.3<0, Term.3, 0) + Term.4
    DBeta[4,]   <-  ifelse(Term.3>0, Term.3, 0)

    if(components=="two"){

      Res        <-  c(sum(DBeta[1:2,]), sum(DBeta[3:4,]))
      names(Res) <-  c("Extinction component", "Colonization component")

    } else if(components=="four"){

      Res        <-  rowSums(DBeta)
      names(Res) <-  c("Extinction homogenization", "Extinction differentiation",
                       "Colonization homogenization", "Colonization differentiation")

    } else if(components=="sp"){

      Res           <-  DBeta
      rownames(Res) <-  c("Extinction homogenization", "Extinction differentiation",
                          "Colonization homogenization", "Colonization differentiation")
      colnames(Res) <-  colnames(d1)

    }
  # End of index=="whittaker"

  } else if(index=="baselga"){

    not.zero <-  function(x) sapply(1:length(x), function(i) all.equal(as.numeric(x[i]), 0)) != TRUE
    A        <-  sum(d1)
    M        <-  apply(d1, 2, max) |> sum()
    Y        <-  N*M / (N-1) / sum(d2)
    DBeta    <-  matrix(nrow=4, ncol=S)

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

      Res         <-  c(sum(DBeta[1:2,]), sum(DBeta[3:4,]))
      names(Res)  <-  c("Subtractive component", "Additive component")

    } else if(components=="four"){

      Res         <-  rowSums(DBeta)
      names(Res)  <-  c("Subtractive homogenization", "Subtractive differentiation",
                        "Additive homogenization", "Additive differentiation")

    } else if(components=="sp"){
      Res           <-  DBeta
      rownames(Res) <-  c("Subtractive homogenization", "Subtractive differentiation",
                          "Additive homogenization", "Additive differentiation")
      colnames(Res) <-  colnames(d1)
    }

  } # End of index=="baselga"

  return(Res)
}
