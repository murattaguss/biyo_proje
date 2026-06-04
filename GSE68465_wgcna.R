# WGCNA - genleri birlikte ifade edilen modullere ayiriyoruz
library(WGCNA)
library(dplyr)
library(tidyverse)

enableWGCNAThreads() # yoksa cok yavas calisiyo

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")

if (!all(rownames(GSE68465) == rownames(df_label))) {
  df_label <- df_label[match(rownames(GSE68465), rownames(df_label)), , drop = FALSE]
}

# 20bin+ gen cok fazla, en degisken 4000 geni aliyoruz
cat("top 4000 gen seciliyor...\n")
gene_vars <- apply(GSE68465, 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:4000])
df <- GSE68465[, top_genes]
df[] <- sapply(df, as.numeric)

# soft threshold secimi - scale-free topology icin en uygun gucu ariyoruz
cat("soft threshold seciliyor...\n")
powers <- c(1:20)
sft <- pickSoftThreshold(df, powerVector = powers, verbose = 5)

png("plots/soft_th.png", width=12, height=6, units="in", res=300)
par(mfrow = c(1, 2))
# soldaki: R^2 ne kadar yuksekse o kadar iyi, 0.8 ustu olsun istiyoruz
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Power", ylab = "Scale Free Topology Model Fit", type = "n",
     main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red", cex = 1)
abline(h = 0.8, col = "red", lty = 2)
# sagdaki: connectivity, guc arttikca duser
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# 0.8'i gecen en kucuk gucu seciyoruz, gecen yoksa 6 aliyoruz
fit_indices <- sft$fitIndices
suitable_powers <- fit_indices$Power[fit_indices$SFT.R.sq > 0.8]
if (length(suitable_powers) > 0) {
  softPower <- min(suitable_powers)
} else {
  softPower <- 6
}
cat("secilen guc:", softPower, "\n")

# adjacency ve TOM matrisleri
cat("adjacency hesaplaniyor...\n")
adjacency <- adjacency(df, power = softPower)

cat("TOM hesaplaniyor...\n")
TOM <- TOMsimilarity(adjacency)
dimnames(TOM) <- list(colnames(df), colnames(df))
dissTOM <- 1 - TOM

# genleri kume agacina gore modullere ayiriyoruz
cat("kumeleme yapiliyor...\n")
geneTree <- hclust(as.dist(dissTOM), method = "average")

minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, cutHeight = 0.99, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
moduleColors <- labels2colors(dynamicMods)
cat("modul boyutlari:\n")
print(table(moduleColors))

png("plots/gene_dendrogram_module.png", width=12, height=6, units="in", res=300)
plotDendroAndColors(geneTree, moduleColors, "Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# benzer modulleri birlestiriyoruz
cat("eigengene'ler hesaplaniyor...\n")
MEList <- moduleEigengenes(df, colors = moduleColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

png("plots/METree.png", width=12, height=6, units="in", res=300)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.25, col = "red")
dev.off()

cat("benzer moduller birlestiriliyor...\n")
merge <- mergeCloseModules(df, moduleColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

png("plots/METree_merged.png", width=12, height=6, units="in", res=300)
plotDendroAndColors(geneTree, cbind(moduleColors, mergedColors),
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

moduleColors <- mergedColors
MEs <- mergedMEs

# her modul icin top 5 hub gen - kME degerine gore seciyoruz
# kME = genin moduldeki eigengene ile korelasyonu, yuksekse hub gen
cat("hub genler seciliyor...\n")
geneModuleMembership <- as.data.frame(cor(df, MEs, use = "p"))
colnames(geneModuleMembership) <- sub("^ME", "", colnames(geneModuleMembership))

all_modules <- unique(moduleColors)
all_modules <- all_modules[all_modules != "grey"] # gri = atanamamis genler, isimize yaramaz

hub_list <- list()
for (mod in all_modules) {
  genes_in_mod <- colnames(df)[moduleColors == mod]
  kme_vals <- geneModuleMembership[genes_in_mod, mod]
  
  mod_df <- data.frame(
    Gene = genes_in_mod,
    kME = kme_vals,
    abs_kME = abs(kme_vals),
    Module = mod,
    stringsAsFactors = FALSE
  )
  
  mod_df_sorted <- mod_df[order(-mod_df$abs_kME), ]
  top_5 <- head(mod_df_sorted, 5)
  hub_list[[mod]] <- top_5
}

hub_genes_df <- do.call(rbind, hub_list)
rownames(hub_genes_df) <- NULL
write.csv(hub_genes_df, file = "results/hubsInEachModule.csv", row.names = FALSE)

# modul-trait iliskisi heatmap
cat("modul-trait korelasyonlari...\n")
trait <- data.frame(
  vital_status = df_label$event,
  time = df_label$time,
  age = df_label$age,
  stage_late = ifelse(df_label$binary == "Late", 1, 0)
)
rownames(trait) <- rownames(df_label)
trait <- trait[match(rownames(df), rownames(trait)), ]

corMat <- cor(MEs, trait, use = "p")
pMat <- corPvalueStudent(corMat, nrow(df))

png("plots/module-trait2.png", width=12, height=6, units="in", res=300)
labeledHeatmap(
  Matrix = corMat,
  xLabels = colnames(trait),
  yLabels = colnames(MEs),
  textMatrix = round(corMat, 2),
  colors = blueWhiteRed(50),
  main = "Module-Trait Relationships"
)
dev.off()

# cytoscape icin dosya aktarimi - en korele modulu seciyoruz
cat("cytoscape export...\n")
vital_cor <- corMat[, "vital_status"]
non_grey_MEs <- names(vital_cor)[names(vital_cor) != "MEgrey"]

if (length(non_grey_MEs) > 0) {
  most_correlated_ME <- non_grey_MEs[which.max(abs(vital_cor[non_grey_MEs]))]
  export_module <- sub("^ME", "", most_correlated_ME)
} else {
  export_module <- "blue"
}
cat("export edilecek modul:", export_module, "\n")

moduleGenes <- moduleColors == export_module
mod_gene_names <- colnames(df)[moduleGenes]
modTOM <- TOM[moduleGenes, moduleGenes]
dimnames(modTOM) <- list(mod_gene_names, mod_gene_names)

exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste0("results/edges_", export_module, ".txt"),
  nodeFile = paste0("results/nodes_", export_module, ".txt"),
  weighted = TRUE,
  threshold = 0.05,
  nodeNames = mod_gene_names,
  altNodeNames = mod_gene_names,
  nodeAttr = moduleColors[moduleGenes]
)

# modul gen listelerini kaydediyoruz
cat("modul uyelikleri kaydediliyor...\n")
df_modules <- data.frame(
  MODULE = moduleColors,
  GENE = colnames(df),
  stringsAsFactors = FALSE
)
write.table(df_modules, file = "results/WGCNAmodules.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

split_modules <- split(df_modules, df_modules$MODULE)
lapply(names(split_modules), function(module) {
  genes_df <- unique(split_modules[[module]])
  fileName <- paste0("results/modules-", module, ",n=", nrow(genes_df), ".txt")
  write.table(genes_df, file = fileName, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
})

cat("WGCNA bitti!\n")
