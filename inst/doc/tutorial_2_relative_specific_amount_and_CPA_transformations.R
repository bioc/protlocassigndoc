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

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
library(protlocassign)

## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------
data(protNSA_test)
data(totProtAT5)
protNSA <- protNSA_test 
str(protNSA)
totProt <- totProtAT5
round(totProt, digits=4)

## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------
protRSA  <- RSAfromNSA(NSA=protNSA[,1:9],
                         NstartMaterialFractions=6, totProt=totProt)
dim(protRSA)
str(protRSA)

## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------
protRSA <- data.frame(protRSA, protNSA[,10:11])
#note data frame is being overwritten
dim(protRSA)
str(protRSA)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
data(markerListJadot)
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA,
                          markerList=markerListJadot, numDataCols=9)
round(refLocationProfilesNSA, digits=4)


## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------
refLocationProfilesRSA <- RSAfromNSA(NSA=refLocationProfilesNSA, NstartMaterialFractions=6,
       totProt=totProt)
round(refLocationProfilesRSA, digits=4)


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesRSA_2 <- locationProfileSetup(profile=protRSA,
                          markerList=markerListJadot, numDataCols=9)
round(refLocationProfilesRSA_2, digits=4)
# we use the `as.matrix` function for display purposes in the tutorial
as.matrix(all.equal(refLocationProfilesRSA, refLocationProfilesRSA_2, 
                    precision=0, countEQ=TRUE))

## ---- echo=TRUE, eval=TRUE, fig.width = 7, fig.height = 11------------------------------------------------------------------------------------------
loc.list <- rownames(refLocationProfilesRSA)
n.loc <- length(loc.list)
par(mfrow=c(4,2))
for (i in 1:n.loc) {
  markerProfilePlot(refLoc=loc.list[i], profile=protRSA,
                     markerList=markerListJadot,
                     refLocationProfiles=refLocationProfilesRSA, ylab="RSA")
  }

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protCPAfromRSA <- fitCPA(profile=protRSA,
                        refLocationProfiles=refLocationProfilesRSA, 
                        numDataCols=9)
str(protCPAfromRSA)

## ---- echo=TRUE, eval=TRUE, fig.width=7, fig.height=11----------------------------------------------------------------------------------------------

protPlotfun(protName="TLN1", profile=protRSA, numDataCols=9,
                        refLocationProfiles=refLocationProfilesRSA,
                        assignPropsMat=protCPAfromRSA,
                        yAxisLabel="Relative Specific Amount")

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

