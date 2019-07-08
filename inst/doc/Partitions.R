## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- eval=TRUE----------------------------------------------------------
library(partitions)
data(Italy)

## ---- eval=TRUE----------------------------------------------------------
tail(Italy, rownames=TRUE)
tail(Id_Italy)

## ---- eval=TRUE----------------------------------------------------------
PosNB <- SampleCluster(data=Italy, Prior="ESCNB", burn=5, nsamples=10, truncds=0.01)

## ---- eval=TRUE----------------------------------------------------------
PosNB$Z[1:5,1:10]
head(PosNB$Params)

## ---- eval=FALSE---------------------------------------------------------
#  PosNB <- SampleCluster(data=Italy, Prior="ESCNB", burn=5, nsamples=10, truncds=0.01, samds = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  PosDP <- SampleCluster(Italy, "DP", 5, 10, 0.01)
#  PosPY <- SampleCluster(Italy, "PY", 5, 10, 0.01)

## ---- eval=TRUE----------------------------------------------------------
PosD <- SampleClusterESCD(Italy, burn=5, nsamples=10, truncds=0.01)
PosD$Z[1:5,1:10]
head(PosD$Params)

## ---- eval=TRUE----------------------------------------------------------
calcError(PosD$Z, Id_Italy$id)

