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
data(protNSA_AT5tmtMS2)
dim(protNSA_AT5tmtMS2)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
data(protNSA_test)
dim(protNSA_test)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
protNSA <- protNSA_test

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
round(head(protNSA), digits=2)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
data(markerListJadot)
dim(markerListJadot)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
markerListJadot

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA, markerList=markerListJadot, numDataCols=9)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
round(refLocationProfilesNSA, digits=4)

## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------
rownames(refLocationProfilesNSA)

## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------

markerProfilePlot(refLoc="PM", profile=protNSA, markerList=markerListJadot,
                     refLocationProfiles=refLocationProfilesNSA, ylab="NSA")

## ---- echo=TRUE, eval=TRUE, fig.show='hold'---------------------------------------------------------------------------------------------------------

markerProfilePlot(refLoc="PM", profile=protNSA, 
                  markerList=markerListJadot,
                  refLocationProfiles=refLocationProfilesNSA, ylab="NSA",
                  refProtPlot=2)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 5, fig.height = 7--------------------------------------------------------------------------

par(mfrow=c(3,2)) # this will be new default layout for subsequent plots.  
# Will need to reset par(mfrow=c(1,1)) for single graph layouts

for (j in 1:5) {
  markerProfilePlot(refLoc="PM", profile=protNSA,
                    markerList=markerListJadot,
                    refLocationProfiles=refLocationProfilesNSA, 
                    ylab="NSA", refProtPlot=j)
  }

## ---- echo=TRUE, eval=TRUE, fig.width = 7, fig.height = 10------------------------------------------------------------------------------------------
loc.list <- rownames(refLocationProfilesNSA)
n.loc <- length(loc.list)
par(mfrow=c(4,2))
for (i in 1:n.loc) {
  markerProfilePlot(refLoc=loc.list[i], profile=protNSA,
                     markerList=markerListJadot,
                     refLocationProfiles=refLocationProfilesNSA, ylab="NSA")
  }

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------
protCPAfromNSA <- fitCPA(profile=protNSA,
                          refLocationProfiles=refLocationProfilesNSA, 
                          numDataCols=9)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
round(tail(protCPAfromNSA), digits=2)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protIndex("TLN1", profile=protNSA)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protIndex("TLN", profile=protNSA)

## ---- echo=TRUE, eval=TRUE, fig.width=7, fig.height=9-----------------------------------------------------------------------------------------------

protPlotfun(protName="TLN1", profile=protNSA, 
numDataCols=9, n.compartments=8,
  refLocationProfiles=refLocationProfilesNSA, 
  assignPropsMat=protCPAfromNSA, 
  yAxisLabel="Normalized Specific Amount")

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

