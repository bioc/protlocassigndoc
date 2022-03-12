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
data(protNSA_test)
dim(protNSA_test)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
data(totProtAT5)
protAcup_test_data <- AcupFromNSA(protNSA_test[,1:6],
                                       NstartMaterialFractions=6, 
                                       totProt=totProtAT5[1:6])

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
# take the first six fractions:
protAcupSixTemp <- protAcup_test_data[,1:6]  
# add L1 and L2 together, and put into L1; then re-name L1 as L
# make L1 the sum of L1 and L2:
protAcupFiveTemp <- within(protAcupSixTemp, L1 <- L1+L2) 
# L is now the name of third columm:
names(protAcupFiveTemp)[3] <- "L"  
# drop L2 column (column 4)
protAcup <- protAcupFiveTemp[, -4]  
head(round(protAcup,3))

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
round(totProtAT5, 3)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
totProtTemp <- totProtAT5
totProtTemp[3] <- totProtAT5[3] + totProtAT5[4]
names(totProtTemp)[3] <- "L"
# remove L2 (component 4) and Nyc values (columns 7, 8, 9)
totProt <- totProtTemp[c(-4, -7, -8, -9)]   
round(totProt, 3)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
protNSA<- NSAfromAcup(Acup=protAcup, NstartMaterialFractions=5, 
                      totProt=totProt)
head(round(protNSA,3))

protRSA <- RSAfromNSA(NSA=protNSA,NstartMaterialFractions = 5, 
                      totProt = totProt)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
data(markerListJadot)
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA, 
                            markerList=markerListJadot, numDataCols=5)
refLocationProfilesAcup <- AcupFromNSA(refLocationProfilesNSA, 
                            NstartMaterialFractions=5, totProt=totProt)
refLocationProfilesRSA <- RSAfromNSA(refLocationProfilesNSA, 
                            NstartMaterialFractions=5, totProt=totProt)

round(refLocationProfilesNSA, 3)
round(refLocationProfilesRSA, 3)
round(refLocationProfilesAcup, 3)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7.5, fig.height = 10-----------------------------------------------------------------------

loc.list <- rownames(refLocationProfilesRSA)
n.loc <- length(loc.list)
par(mfrow=c(4,2))
for (i in 1:n.loc) {
  markerProfilePlot(refLoc=loc.list[i], profile=protRSA,
                    markerList=markerListJadot,
                    refLocationProfiles=refLocationProfilesRSA, 
                    ylab="RSA")
}

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------

i=1  #Cyto
j=2  #ER
mixProtiProtj12Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj12RSA <- RSAfromAcup(Acup=mixProtiProtj12Acup, 
                            NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA12 <- fitCPA(profile=mixProtiProtj12RSA,
                            refLocationProfiles=refLocationProfilesRSA,
                            numDataCols=5)

i=1  #Cyto
j=3  #Golgi
mixProtiProtj13Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj13RSA <- RSAfromAcup(Acup=mixProtiProtj13Acup, 
                                  NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA13 <- fitCPA(profile=mixProtiProtj13RSA,
                                 refLocationProfiles=refLocationProfilesRSA,
                                 numDataCols=5)
i=1  #Cyto
j=4  #Lyso
mixProtiProtj14Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj14RSA <- RSAfromAcup(Acup=mixProtiProtj14Acup, 
                            NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA14 <- fitCPA(profile=mixProtiProtj14RSA,
                                refLocationProfiles=refLocationProfilesRSA,
                                numDataCols=5)

i=1  #Cyto
j=5  #Mito
mixProtiProtj15Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj15RSA <- RSAfromAcup(Acup=mixProtiProtj15Acup, 
                            NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA15 <- fitCPA(profile=mixProtiProtj15RSA,
                                refLocationProfiles=refLocationProfilesRSA,
                                numDataCols=5)

i=1  #Cyto
j=6  #Nuc
mixProtiProtj16Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj16RSA <- RSAfromAcup(Acup=mixProtiProtj16Acup, 
                            NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA16 <- fitCPA(profile=mixProtiProtj16RSA,
                                refLocationProfiles=refLocationProfilesRSA,
                                numDataCols=5)
i=1  #Cyto
j=7  #Perox
mixProtiProtj17Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj17RSA <- RSAfromAcup(Acup=mixProtiProtj17Acup, 
                            NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA17 <- fitCPA(profile=mixProtiProtj17RSA,
                                refLocationProfiles=refLocationProfilesRSA,
                                numDataCols=5)

i=1  #Cyto
j=8  #PM
mixProtiProtj18Acup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)
mixProtiProtj18RSA <- RSAfromAcup(Acup=mixProtiProtj18Acup, 
                            NstartMaterialFractions=5, totProt=totProt)
mixProtiProtjCPAfromRSA18 <- fitCPA(profile=mixProtiProtj18RSA,
                             refLocationProfiles=refLocationProfilesRSA,
                            numDataCols=5)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7.5, fig.height = 10-----------------------------------------------------------------------

#library(pracma)   # for trapz function
par(mfrow=c(2,4))

mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA12,
            NstartMaterialFractions=5, Loc1=1, Loc2=2,
            errorReturn = TRUE, subTitle="RSA")


mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA13,
            NstartMaterialFractions=5, Loc1=1, Loc2=3,
            errorReturn = TRUE, subTitle="RSA")

mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA14,
            NstartMaterialFractions=5, Loc1=1, Loc2=4,
            errorReturn = TRUE, subTitle="RSA")

mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA15,
            NstartMaterialFractions=5, Loc1=1, Loc2=5,
            errorReturn = TRUE, subTitle="RSA")

mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA16,
            NstartMaterialFractions=5, Loc1=1, Loc2=6,
            errorReturn = TRUE, subTitle="RSA")

mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA17,
            NstartMaterialFractions=5, Loc1=1, Loc2=7,
            errorReturn = TRUE, subTitle="RSA")


mixturePlot(mixProtiProtjCPA=mixProtiProtjCPAfromRSA18,
            NstartMaterialFractions=5, Loc1=1, Loc2=8,
            errorReturn = TRUE, subTitle="RSA")


## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7.5, fig.height = 10-----------------------------------------------------------------------
#install.packages(c("plot.matrix", "viridis", "grid", "gridExtra"))

par(mfrow=c(1,1))
errorMatAll <- mixtureHeatMap(Acup=refLocationProfilesAcup, 
                        totProt=totProt, NstartMaterialFractions=5)


## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 2--------------------------------------------------------------------------
errorMatAll
#library(plot.matrix)
#library("viridis")
op <- par(mar=c(4,4,1,1))  # save default parameters in variable op

col <- rev(viridis::magma(16))
library(plot.matrix)
plot(errorMatAll, col=col, breaks=seq(0, 15, 1), key=NULL, main="",
     axis.col=NULL, axis.row=NULL, 
     xlab="RSA                        NSA                       Acup", 
     ylab="Log2 Linear", digits=2, cex=2, cex.lab=1.1)
par(op) # restore default parameters


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------

data.frame(rownames(mixProtiProtj14RSA))
mixCyto50Lyso50RSA <- mixProtiProtj14RSA[6,]


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
mixCyto50Lyso50CPAfromRSAforced <- fitCPA(profile= mixCyto50Lyso50RSA,
                                refLocationProfiles=refLocationProfilesRSA,
                                numDataCols=5, ind.vary=c(1,4), minVal=TRUE)
round(mixCyto50Lyso50CPAfromRSAforced, digits=4)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
mixCyto50Lyso50CPAfromRSApairs <- fCPAsubsets(profile= mixCyto50Lyso50RSA,
                                  refLocationProfiles=refLocationProfilesRSA,
                                  numDataCols=5,nCPAcomparts=2)
round(mixCyto50Lyso50CPAfromRSApairs, digits=4)

## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
mixCyto50Lyso50CPAfromRSAallSubsets<- NULL
for (kk in 1:8) {
   temp  <- fCPAsubsets(profile=mixCyto50Lyso50RSA,
                              refLocationProfiles=refLocationProfilesRSA,
                              numDataCols=5, nCPAcomparts=kk)
   mixCyto50Lyso50CPAfromRSAallSubsets <- 
                       rbind(mixCyto50Lyso50CPAfromRSAallSubsets, temp)
}
# sort dataframe by the variable “value” by creating an intermediate 
#   list of indices "ord_fit".
ord_fit <- order(mixCyto50Lyso50CPAfromRSAallSubsets$value)
mixCyto50Lyso50CPAfromRSAallSubsetsOrd <- 
         mixCyto50Lyso50CPAfromRSAallSubsets[ord_fit,]
round(head(mixCyto50Lyso50CPAfromRSAallSubsetsOrd), 3)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 5--------------------------------------------------------------------------
modelVals <- 1:255
par(mfrow=c(1,1))  # reset the plot space
plot(mixCyto50Lyso50CPAfromRSAallSubsetsOrd$value ~ modelVals, 
     xlab="Model id", ylab="Goodness of fit value", las=1)

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7, fig.height = 5--------------------------------------------------------------------------

par(mfrow=c(1,1))  # reset the plot space
plot(mixCyto50Lyso50CPAfromRSAallSubsetsOrd$Lyso ~ modelVals, xlab="Model id", 
     ylab="CPA proortion", type="n", ylim=c(0, 1.1), las=1)

points(mixCyto50Lyso50CPAfromRSAallSubsetsOrd$Lyso ~ modelVals, pch=4, 
       col="darkgreen", cex=1.5)
points(mixCyto50Lyso50CPAfromRSAallSubsetsOrd$Cyto ~ modelVals, pch=1, col="red")
legend(x="topleft", legend=c("Cyto proportion", "Lyso proportion"),
       pch=c( 1, 4), col=c("red", "darkgreen"))

## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7.5, fig.height = 10-----------------------------------------------------------------------


errorAllRSAlinear <- mixturePlotPanel(refLocationProfilesAcup=
                                refLocationProfilesAcup,
                                totProt=totProt, NstartMaterialFractions=5, 
                                errorReturn = TRUE,
                                fitType="RSA", log2Transf=FALSE)


## ---- echo=TRUE, eval=TRUE, fig.show='hold', fig.width = 7.5, fig.height = 10-----------------------------------------------------------------------

errorAllRSAlog2 <- mixturePlotPanel(refLocationProfilesAcup=
                                    refLocationProfilesAcup,
                                    totProt=totProt, 
                                    NstartMaterialFractions=5, 
                                    errorReturn = TRUE,
                                    fitType="RSA", log2Transf=TRUE)


errorAllNSAlinear <- mixturePlotPanel(refLocationProfilesAcup=
                                    refLocationProfilesAcup,
                                    totProt=totProt, 
                                    NstartMaterialFractions=5, 
                                    errorReturn = TRUE,
                                    fitType="NSA", log2Transf=FALSE)

errorAllNSAlog2 <- mixturePlotPanel(refLocationProfilesAcup=
                                    refLocationProfilesAcup,
                                    totProt=totProt, 
                                    NstartMaterialFractions=5, 
                                    errorReturn = TRUE,
                                    fitType="NSA", log2Transf=TRUE)

errorAllAcupLinear <- mixturePlotPanel(refLocationProfilesAcup=
                                    refLocationProfilesAcup,
                                    totProt=totProt, 
                                    NstartMaterialFractions=5, 
                                    errorReturn = TRUE,
                                    fitType="Acup", log2Transf=FALSE)

errorAllAcupLog2 <- mixturePlotPanel(refLocationProfilesAcup=
                                    refLocationProfilesAcup,
                                    totProt=totProt, 
                                    NstartMaterialFractions=5, 
                                    errorReturn = TRUE,
                                    fitType="Acup", log2Transf=TRUE)



## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------
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

