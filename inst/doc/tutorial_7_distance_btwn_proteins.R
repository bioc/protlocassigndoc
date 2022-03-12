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
data(protNSA_AT5tmtMS2)
data(totProtAT5)
protNSA <- protNSA_AT5tmtMS2
totProt <- totProtAT5
distUseNSA <- dist(protNSA[,1:9], method="euclidean")

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protsUse <- rownames(protNSA)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
nearestProts(protName="CTSD", n.nearest=10,  distProts=distUseNSA, protNames=protsUse,
             profile=protNSA)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protProfileLevelsRSA <- RSAfromNSA(NSA=protNSA[,1:9],
                                 NstartMaterialFractions=6, totProt=totProt)
distUseRSA <- dist(protProfileLevelsRSA, method="euclidean")
nearestProts(protName="CTSD", n.nearest=10,  distProts=distUseRSA, protNames=protsUse,
             profile=protProfileLevelsRSA)


## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
distUseNSAmatrix <- as.matrix(distUseNSA)
distUseNSAmatrix[1:5,1:5]

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

