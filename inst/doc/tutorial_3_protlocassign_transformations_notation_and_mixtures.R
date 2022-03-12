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
data(markerListJadot)
protNSA <- protNSA_test
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA, 
                           markerList=markerListJadot, numDataCols=9)
round(refLocationProfilesNSA, digits=3)


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
data(totProtAT5)
totProt <- totProtAT5
round(totProt, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
sum(totProt[1:6])

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------

refLocationProfilesAcup <- AcupFromNSA(NSA=refLocationProfilesNSA, NstartMaterialFractions=6, 
                               totProt=totProt)
round(refLocationProfilesAcup, digits=4)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesRSA <- RSAfromAcup(refLocationProfilesAcup, NstartMaterialFractions=6, totProt=totProt)
round(refLocationProfilesRSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesRSA_2 <- RSAfromNSA(NSA=refLocationProfilesNSA,
                                NstartMaterialFractions=6, totProt=totProt)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
identical(refLocationProfilesRSA, refLocationProfilesRSA_2)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesNSA_2 <- NSAfromRSA(refLocationProfilesRSA_2)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
as.matrix(all.equal(refLocationProfilesNSA_2, refLocationProfilesNSA, 
                    tolerance=0, countEQ=TRUE))

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------

refLocationProfilesAcup <- AcupFromNSA(NSA=refLocationProfilesNSA, 
                                       NstartMaterialFractions=6, 
                                       totProt=totProt)
data.frame(rownames(refLocationProfilesAcup))
mixCytoLysoAcup <- proteinMix(AcupRef=refLocationProfilesAcup, 
                              increment.prop=0.1,
                              Loc1=1, Loc2=4)
# Note that the default value of increment.prop=0.1. 
# This does not need to be explicitly 
#  specified unless a different increment is desired.

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
round(mixCytoLysoAcup, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixCytoLysoRSA <- RSAfromAcup(Acup=mixCytoLysoAcup, 
                              NstartMaterialFractions=6, totProt=totProt)

round(mixCytoLysoRSA, digits=3)


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixCytoLysoCPAfromRSA <- fitCPA(profile=mixCytoLysoRSA,
                             refLocationProfiles=refLocationProfilesRSA,
                            numDataCols=9)
round(mixCytoLysoCPAfromRSA, digits=3)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
library(pracma)
par(mfrow=c(1,1))  # reset window for a single plot
# The argument increment.prop needs to match the value used  
#   in creating the mixture using proteinMix. This does
#   not need to be specified if using the value of 0.1
mixturePlot(mixProtiProtjCPA=mixCytoLysoCPAfromRSA, 
            NstartMaterialFractions=6, Loc1=1, Loc2=4,
            increment.prop=0.1, xaxisLab=TRUE, yaxisLab=TRUE)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 10-------------------------------------------------------------------------
par(mfrow=c(4,2))
i <- 1
mixErrorMat <- NULL
for (j in 2:8) {   
   # Create the mixture of Cyto (i = 1) with compartment j
   mixProtiProtjAcup <- proteinMix(Acup=refLocationProfilesAcup, 
                                   Loc1=i, Loc2=j)

   # Tranform the mixtures to relative specific amounts
   mixProtiProtjRSA <- RSAfromAcup(Acup=mixProtiProtjAcup, 
                        NstartMaterialFractions=6, totProt=totProt)
    
   # Find the constrained proportional assigments (CPA)                 
   mixProtiProtjCPAfromRSA <- fitCPA(profile=mixProtiProtjRSA,
                    refLocationProfiles=refLocationProfilesRSA, 
                    numDataCols=9)
   
   # Plot the results, including the area-based error estimate, 
   #    and collect the area-based errors (errorReturn=TRUE)                         
    mixResult <- mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA, 
                             NstartMaterialFractions=6, Loc1=i, Loc2=j, 
                             increment.prop=0.1, errorReturn = TRUE)
    mixErrorMat <- rbind(mixErrorMat, mixResult)         
}
mixErrorAllCytoRSA <- mixErrorMat

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
mixErrorAllCytoRSA

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

