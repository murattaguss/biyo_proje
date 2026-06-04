
# projede lazim olan butun paketleri buraya topladim, her seferinde teker teker kurmakla ugrasmasin diye
# BUG COZUCU CRAN SABITLEME
options(repos = c(CRAN = "https://cloud.r-project.org"))

# bioconductor paketleri icin once bu lazim
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# raw veri isleme icin gerekli
BiocManager::install("GEOquery", update = FALSE)
BiocManager::install("affy", update = FALSE)
BiocManager::install("hgu133a.db", update = FALSE)

# limma - DGE analizi icin, hoca bunu kullanin demisti
BiocManager::install("limma", update = FALSE)
library(limma)

BiocManager::install("Biobase", update = FALSE)
library(Biobase)

# wgcna icin bunlar lazimmis, eksik olunca hata veriyo
BiocManager::install("impute", update = FALSE)
library(impute)

BiocManager::install("preprocessCore", update = FALSE)
library(preprocessCore)

install.packages("WGCNA")
library(WGCNA)

# GO analizi icin
BiocManager::install("clusterProfiler", update = FALSE)
library(clusterProfiler)

BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)

BiocManager::install("AnnotationDbi", update = FALSE)
BiocManager::install("GO.db", update = FALSE)
BiocManager::install("enrichplot", update = FALSE)

# genel veri isleme ve grafik paketleri
install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)

install.packages("ggplot2")
library(ggplot2)
install.packages("reshape2")
library(reshape2)

# boxplotlara p degeri yazdirmak icin
install.packages("ggpubr")
library(ggpubr)

# survival analizi
install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)

# makine ogrenmesi paketleri
install.packages("randomForest")
library(randomForest)

install.packages("e1071") # svm bunun icinde
library(e1071)

install.packages("glmnet") # lasso icin
library(glmnet)

install.packages("pROC") # roc egrisi cizdirmek icin
library(pROC)

install.packages("gridExtra") # confusion matrix grid icin
library(gridExtra)
