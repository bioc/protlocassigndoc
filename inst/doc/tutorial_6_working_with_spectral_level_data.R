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
library(outliers)
library(lme4)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
library(protlocassign)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
data(spectraNSA_test)
spectraNSA <- spectraNSA_test

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
eps <- 0.029885209
library(BiocParallel)
flagSpectraBox <- outlierFind(protClass=spectraNSA, 
              outlierLevel="peptide", numRefCols=5, numDataCols=9, 
              outlierMeth="boxplot", range=3, eps=eps, 
              randomError=TRUE, setSeed=17356, cpus=4)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
flagSpectraScore <- outlierFind(protClass=spectraNSA, 
                outlierLevel="peptide", numRefCols=5, numDataCols=9, 
                outlierMeth="scores", proba=0.99,
                eps=eps, randomError=TRUE, setSeed=27681, cpus=4)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
str(flagSpectraBox, strict.width="cut", width=65)
str(flagSpectraScore, strict.width="cut", width=65)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
length(flagSpectraBox$outlier.num.spectra)    # total number of spectra
sum(flagSpectraBox$outlier.num.spectra == 0)  # number of acceptable spectra
sum(flagSpectraBox$outlier.num.spectra != 0)  # number of outlier spectra

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
length(flagSpectraScore$outlier.num.spectra)    # total number of spectra
sum(flagSpectraScore$outlier.num.spectra == 0)  # number of acceptable spectra
sum(flagSpectraScore$outlier.num.spectra != 0)  # number of outlier spectra

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
flagSpectraBoxRange110ranErrF <- outlierFind(protClass=spectraNSA, 
              outlierLevel="peptide", numRefCols=5, numDataCols=9, 
              outlierMeth="boxplot", range=110, eps=eps, 
              randomError=FALSE, setSeed=NULL, cpus=4)
sum(flagSpectraBoxRange110ranErrF$outlier.num.spectra != 0)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
flagSpectraBoxRange110ranErrT <- outlierFind(protClass=spectraNSA, 
              outlierLevel="peptide", numRefCols=5, numDataCols=9, 
              outlierMeth="boxplot", range=110, eps=eps, 
              randomError=TRUE, setSeed=652908, cpus=4)
sum(flagSpectraBoxRange110ranErrT$outlier.num.spectra != 0)

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
pepProfiles <- profileSummarize(protsCombineCnew=flagSpectraBox,
                numRefCols=6, numDataCols=9, refColsKeep=c(1,2,4),eps=eps,
                GroupBy="peptideId", outlierExclude="spectra", cpus=4)
str(pepProfiles, strict.width="cut", width=65)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
flagPepsBox <- outlierFind(protClass=pepProfiles,            
                    outlierLevel="protein", 
                    numRefCols=3, numDataCols=9, eps=eps,
                    setSeed=823537, cpus=4)
str(flagPepsBox, strict.width="cut", width=65)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
length(flagPepsBox$outlier.num.peptides)
sum(flagPepsBox$outlier.num.peptides == 0)
sum(flagPepsBox$outlier.num.peptides != 0)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
data.frame(names(flagPepsBox))
flagSpectraPeps <- merge(x=flagSpectraBox, 
                        y=flagPepsBox[,c(4,17)], by="pepId")
str(flagSpectraPeps,  strict.width="cut", width=65)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
length(flagSpectraPeps$outlier.num.peptides)
sum({flagSpectraPeps$outlier.num.spectra == 0} & 
      {flagSpectraPeps$outlier.num.peptides == 0})

## ---- echo=TRUE, message=TRUE, warning=TRUE---------------------------------------------------------------------------------------------------------
protNSA_1 <- 
          profileSummarize(protsCombineCnew=flagSpectraPeps,
          numRefCols=7, numDataCols=9, eps=eps, GroupBy="protId",
          outlierExclude="spectraAndPeptide", cpus=4, singularList=TRUE)
str(protNSA_1,  strict.width="cut", width=65)

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
protNSA_2 <- 
          profileSummarize(protsCombineCnew=flagSpectraPeps,
          numRefCols=7, numDataCols=9, eps=eps, GroupBy="protId",
          outlierExclude="spectra", cpus=4)
str(protNSA_2,  strict.width="cut", width=65)

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
protNSA_3 <- 
          profileSummarize(protsCombineCnew=flagSpectraPeps,
          numRefCols=7, numDataCols=9, eps=eps, GroupBy="protId",
          outlierExclude="none", cpus=4)
str(protNSA_3,  strict.width="cut", width=65)

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
sum((protNSA_3$Nspectra))
sum((protNSA_2$Nspectra))
sum((protNSA_1$Nspectra))

sum((protNSA_3$Npep))
sum((protNSA_2$Npep))
sum((protNSA_1$Npep))

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
sum(is.na(protNSA_1$M))
sum(is.na(protNSA_2$M))
sum(is.na(protNSA_3$M))

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------
protNSA_new <- protNSA_1 
#copy first field in data frame containing protein names to rownnames field
rownames(protNSA_new) <- protNSA_1[,1]  
protNSA_new <- protNSA_new[,-1]  #delete first internal column containing protein names

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------
data(protNSA_test)
all.equal(protNSA_test, protNSA_new, countEQ=TRUE, tolerance=0)

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------
data.frame(names(flagPepsBox))
protPepProfileNSA <- protPepProfile(flagPeps=flagPepsBox,
                                 numRefCols=4, numDataCols=9, 
                                 protProfileData=protNSA_1)
str(protPepProfileNSA, strict.width="cut", width=65)

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------
rownames(protPepProfileNSA) <- protPepProfileNSA$peptide

## ---- echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------
head(rownames(protPepProfileNSA))

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
library(protlocassign)

data(totProtAT5)
protPepNSA <- protPepProfileNSA
str(protPepNSA, strict.width="cut", width=65)
totProt <- totProtAT5
totProt

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
data(markerListJadot)
refLocationProfilesNSA <- locationProfileSetup(profile=protPepNSA[, 4 + (1:9)],
                          markerList=markerListJadot, numDataCols=9)
round(refLocationProfilesNSA, digits=4)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
refLocationProfilesRSA <- RSAfromNSA(NSA=refLocationProfilesNSA,
                              NstartMaterialFractions=6, totProt=totProtAT5)
round(refLocationProfilesRSA, digits=4)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protPepRSA_trimmed <- RSAfromNSA(NSA=protPepNSA[, 4 + (1:9)],
                              NstartMaterialFractions=6, totProt=totProtAT5)
str(protPepRSA_trimmed, strict.width="cut", width=65)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protPepRSA <- data.frame(protPepNSA[, 1:4], protPepRSA_trimmed, protPepNSA[,14:15] )  # add in the ref columns
str(protPepRSA, strict.width="cut", width=65)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protRSA.ind <- {protPepRSA$prot == protPepRSA$peptide}  # protein indicators
protRSA <- protPepRSA[protRSA.ind,]  # these are the data for proteins only
dim(protRSA)
pepRSA <- protPepRSA[!protRSA.ind,] # these are the data for peptides only

data.frame(colnames(protRSA))

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protCPAfromRSA <-  fitCPA(profile=protRSA[, 4+1:9],
                      refLocationProfiles=refLocationProfilesRSA, 
                      numDataCols=9)
str(protCPAfromRSA, strict.width="cut", width=65)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7.5, fig.height = 10-----------------------------------------------------------------------
#windows(width=7.5, height=10)  # open a window 7.5 by 10 inches
protPepPlotfun(protName="TLN1", protProfile=protRSA[,5:15],
               Nspectra=TRUE, pepProfile=pepRSA, numRefCols=4,
               numDataCols=9, n.compartments=8, 
               refLocationProfiles=refLocationProfilesRSA,
               assignPropsMat=protCPAfromRSA, 
               yAxisLabel="Relative Specific Amount")

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
protPepCPAfromRSA <- fitCPA(profile=protPepRSA[,4 + 1:9],
                               refLocationProfiles=refLocationProfilesRSA, numDataCols=9)
str(protPepCPAfromRSA, strict.width="cut", width=65)

## ---- echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------------------------------------------
utils::sessionInfo()

