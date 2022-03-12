## ----setup, include = FALSE-------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  #collapse = TRUE,
  comment = "#>",
  fig.width = 4,
  fig.height = 4,
  message = FALSE,
  warning = FALSE,
  tidy.opts = list(
    keep.blank.line = TRUE,
    width.cutoff = 150
  ),
  options(width = 150),
  eval = TRUE
)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
library(protlocassign)
data(protNSA_test)
data(totProtAT5)
protNSA <- protNSA_test
totProt <- totProtAT5

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
data(markerListJadot)
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA, 
                            markerList=markerListJadot, numDataCols=9)
refLocationProfilesAcup <- AcupFromNSA(refLocationProfilesNSA, 
                            NstartMaterialFractions=6, totProt=totProt)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
i=1
j=4
mixProt1Prot4Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProt1Prot4RSA <- RSAfromAcup(Acup=mixProt1Prot4Acup, 
                            NstartMaterialFractions=6, totProt=totProt)
# increment.prop is not specified so we use the default value of 0.1


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
round(mixProt1Prot4Acup, digits=3)
round(mixProt1Prot4RSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesRSA <- RSAfromNSA(NSA=refLocationProfilesNSA, 
                          NstartMaterialFractions=6,
                          totProt=totProt)
mixProt1Prot4CPAfromRSA <- fitCPA(profile=mixProt1Prot4RSA, 
                      refLocationProfiles=refLocationProfilesRSA, 
                      numDataCols=9)
round(mixProt1Prot4CPAfromRSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixProt1Prot4NSA <- NSAfromRSA(mixProt1Prot4RSA)
round(mixProt1Prot4NSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------

mixProt1Prot4CPAfromNSA <- fitCPA(profile=mixProt1Prot4NSA,
                             refLocationProfiles=refLocationProfilesNSA,
                             numDataCols=9)
round(mixProt1Prot4CPAfromNSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixProt1Prot4CPAfromAcup <- 
                  fitCPA(profile=mixProt1Prot4Acup,
                          refLocationProfiles=refLocationProfilesAcup, 
                         numDataCols=9)
round(mixProt1Prot4CPAfromAcup, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
round(refLocationProfilesAcup, digits=4)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 4--------------------------------------------------------------------------
par(mfrow=c(1,3))
# In the following, the argument increment.prop is not specified so
#  we use the default value of 0.1
mixturePlot(mixProtiProtjCPA=mixProt1Prot4CPAfromRSA, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j, 
            errorReturn = TRUE, subTitle="RSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot4CPAfromNSA, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="NSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot4CPAfromAcup, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j, 
            errorReturn = TRUE, subTitle="Acup")


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
# As discussed earlier, the argument increment.prop=0.1 does not
#  need to be specified because it is the default
#  and is consistent with the previously generated mixtures.
mixtureAreaError(mixProtiProtjCPA=mixProt1Prot4CPAfromRSA, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j)
mixtureAreaError(mixProtiProtjCPA=mixProt1Prot4CPAfromNSA, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j)
mixtureAreaError(mixProtiProtjCPA=mixProt1Prot4CPAfromAcup, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
eps <- 0.001


mixProt1Prot4CPAfromRSAlog2 <- fitCPA(profile=log2(mixProt1Prot4RSA + eps), 
                      refLocationProfiles=log2(refLocationProfilesRSA + eps), 
                      numDataCols=9)
round(mixProt1Prot4CPAfromRSAlog2, digits=3)

mixProt1Prot4CPAfromNSAlog2 <- 
            fitCPA(profile=log2(mixProt1Prot4NSA + eps),
                refLocationProfiles=log2(refLocationProfilesNSA + eps),
                numDataCols=9)
round(mixProt1Prot4CPAfromNSAlog2, digits=3)

mixProt1Prot4CPAfromAcupLog2 <- 
            fitCPA(profile=log2(mixProt1Prot4Acup + eps),
              refLocationProfiles=log2(refLocationProfilesAcup + eps),
              numDataCols=9)
round(mixProt1Prot4CPAfromAcupLog2, digits=3)


## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 4--------------------------------------------------------------------------
par(mfrow=c(1,3))
mixturePlot(mixProtiProtjCPA=mixProt1Prot4CPAfromRSAlog2, 
            NstartMaterialFractions=6,
            Loc1=i, Loc2=j, errorReturn = TRUE,
            subTitle="Log2 RSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot4CPAfromNSAlog2, 
            NstartMaterialFractions=6, 
            Loc1=i, Loc2=j, errorReturn = TRUE,
            subTitle="Log2 NSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot4CPAfromAcupLog2, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="Log2 Acup")

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
i=1  # Cyto
j=6  # Nuc
mixProt1Prot6Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProt1Prot6RSA <- RSAfromAcup(Acup=mixProt1Prot6Acup, 
                          NstartMaterialFractions=6, totProt=totProtAT5)


mixProt1Prot6CPAfromRSA <- fitCPA(profile=mixProt1Prot6RSA, 
            refLocationProfiles=refLocationProfilesRSA, numDataCols=9)
round(mixProt1Prot6CPAfromRSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixProt1Prot6NSA <- NSAfromRSA(mixProt1Prot6RSA)
mixProt1Prot6CPAfromNSA <- 
                    fitCPA(profile=mixProt1Prot6NSA,
                    refLocationProfiles=refLocationProfilesNSA, 
                    numDataCols=9)
round(mixProt1Prot6CPAfromNSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixProt1Prot6CPAfromAcup <- fitCPA(profile=mixProt1Prot6Acup,
                             refLocationProfiles=refLocationProfilesAcup, 
                             numDataCols=9)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
round(mixProt1Prot6CPAfromAcup, digits=3)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 4--------------------------------------------------------------------------
par(mfrow=c(1,3))
mixturePlot(mixProtiProtjCPA=mixProt1Prot6CPAfromRSA, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="RSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot6CPAfromNSA, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="NSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot6CPAfromAcup, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="Acup")

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
eps <- 0.001

mixProt1Prot6CPAfromRSAlog2 <- 
           fitCPA(profile=log2(mixProt1Prot6RSA + eps), 
                refLocationProfiles=log2(refLocationProfilesRSA + eps),
                numDataCols=9)
round(mixProt1Prot6CPAfromRSAlog2, digits=3)

mixProt1Prot6CPAfromNSAlog2 <- 
           fitCPA(profile=log2(mixProt1Prot6NSA + eps),
                refLocationProfiles=log2(refLocationProfilesNSA + eps), 
                numDataCols=9)
round(mixProt1Prot6CPAfromNSAlog2, digits=3)


mixProt1Prot6CPAfromAcupLog2 <- 
           fitCPA(profile=log2(mixProt1Prot6Acup + eps),
               refLocationProfiles=log2(refLocationProfilesAcup + eps), 
               numDataCols=9)
round(mixProt1Prot6CPAfromAcupLog2, digits=3)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 4--------------------------------------------------------------------------
par(mfrow=c(1,3))
mixturePlot(mixProtiProtjCPA=mixProt1Prot6CPAfromRSAlog2, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="Log2RSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot6CPAfromNSAlog2, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="Log2NSA")
mixturePlot(mixProtiProtjCPA=mixProt1Prot6CPAfromAcupLog2, 
            NstartMaterialFractions=6, Loc1=i, Loc2=j,
            errorReturn = TRUE, subTitle="Log2Acup")

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 4--------------------------------------------------------------------------
par(mfrow=c(1,1))
errorMatAll <- mixtureHeatMap(Acup=refLocationProfilesAcup, totProt=totProt)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
errorMatAll 

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 2--------------------------------------------------------------------------
#library(plot.matrix)
#library("viridis")
op <- par(mar=c(4,4,1,1))  # save default parameters in variable op

col <- rev(viridis::magma(12))
plot(errorMatAll, col=col, breaks=seq(0, 12, 1), key=NULL, main="",
     axis.col=NULL, axis.row=NULL, 
     xlab="RSA                        NSA                       Acup", 
     ylab="Log2 Identity", digits=2, cex=2, cex.lab=1.1)
par(op) # restore default parameters


## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 11-------------------------------------------------------------------------
# If the function produces a figure margins error message, 
#    open a 7 by 11 window:
# windows(width=7, height=11)

errorAllRSAlinear <- mixturePlotPanel(refLocationProfilesAcup=
                                        refLocationProfilesAcup, 
              totProt=totProtAT5, NstartMaterialFractions=6, errorReturn = TRUE, 
              fitType="RSA", log2Transf=FALSE)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
errorAllRSAlinear

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 11-------------------------------------------------------------------------
sum(errorAllRSAlinear[,3])

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 11-------------------------------------------------------------------------



fitType <- "RSA"
log2Transf <- TRUE


errorAllRSAlog2 <- mixturePlotPanel(refLocationProfilesAcup=
              refLocationProfilesAcup, totProt=totProtAT5, 
              NstartMaterialFractions=6, errorReturn = TRUE, 
              fitType=fitType, log2Transf=log2Transf, eps=0.001)



fitType <- "NSA"
log2Transf <- FALSE

errorAllNSAlinear <- mixturePlotPanel(refLocationProfilesAcup=
                               refLocationProfilesAcup,
              totProt=totProtAT5, NstartMaterialFractions=6,
              errorReturn = TRUE, 
              fitType=fitType, log2Transf=log2Transf)

fitType <- "NSA"
log2Transf <- TRUE


errorAllNSAlog2 <- mixturePlotPanel(refLocationProfilesAcup=
                   refLocationProfilesAcup, totProt=totProtAT5, 
                   NstartMaterialFractions=6, errorReturn = TRUE, 
                   fitType=fitType, log2Transf=log2Transf, eps=0.001)

fitType <- "Acup"
log2Transf <- FALSE


errorAllAcupLinear <- mixturePlotPanel(refLocationProfilesAcup=
              refLocationProfilesAcup, totProt=totProtAT5, 
              NstartMaterialFractions=6, errorReturn = TRUE, 
              fitType=fitType, log2Transf=log2Transf)

fitType <- "Acup"
log2Transf <- TRUE


errorAllAcupLog2 <- mixturePlotPanel(refLocationProfilesAcup=
              refLocationProfilesAcup, totProt=totProtAT5, 
              NstartMaterialFractions=6, errorReturn = TRUE, 
              fitType=fitType, log2Transf=log2Transf, eps=0.001)


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------

errorAllRSAlinear
sum(errorAllRSAlinear[,3])

errorAllRSAlog2
sum(errorAllRSAlog2[,3])

errorAllNSAlinear
sum(errorAllNSAlinear[,3])

errorAllNSAlog2
sum(errorAllNSAlog2[,3])

errorAllAcupLinear
sum(errorAllAcupLinear[,3])

errorAllAcupLog2
sum(errorAllAcupLog2[,3])

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

