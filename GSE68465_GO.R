# GO analizi - WGCNA'dan cikan moduldeki genlerin ne ise yaradigini ogreniyoruz
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(WGCNA)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("veriler yukleniyor...\n")
df_modules <- read.table("results/WGCNAmodules.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")

if (!all(rownames(GSE68465) == rownames(df_label))) {
  df_label <- df_label[match(rownames(GSE68465), rownames(df_label)), , drop = FALSE]
}

df_wgcna <- GSE68465[, df_modules$GENE]
df_wgcna[] <- sapply(df_wgcna, as.numeric)

# vital status ile en korele modulu buluyoruz
cat("hedef modul belirleniyor...\n")
MEList <- moduleEigengenes(df_wgcna, colors = df_modules$MODULE)
MEs <- MEList$eigengenes

trait <- df_label$event
cor_MEs <- cor(MEs, trait, use = "p")

# gri disinda en yuksek korelasyonlu olan
non_grey_MEs <- rownames(cor_MEs)[rownames(cor_MEs) != "MEgrey"]
most_correlated_ME <- non_grey_MEs[which.max(abs(cor_MEs[non_grey_MEs, 1]))]
export_module <- sub("^ME", "", most_correlated_ME)
cat("hedef modul:", export_module, "\n")

target_genes <- df_modules$GENE[df_modules$MODULE == export_module]
cat("gen sayisi:", length(target_genes), "\n")

# enrichGO ile biyolojik surecler (BP) icin zenginlestirme yapiyoruz
cat("GO analizi basliyor...\n")
ego <- enrichGO(
  gene          = target_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  cat("anlamli terimler bulundu, grafikler ciziliyor...\n")
  
  png("plots/GO_barplot.png", width=10, height=8, units="in", res=300)
  p_bar <- barplot(ego, showCategory = 20, title = paste("GO BP Enrichment -", export_module, "Module"))
  print(p_bar)
  dev.off()
  
  png("plots/GO_dotplot.png", width=10, height=8, units="in", res=300)
  p_dot <- dotplot(ego, showCategory = 20, title = paste("GO BP Enrichment -", export_module, "Module"))
  print(p_dot)
  dev.off()
  
  write.csv(as.data.frame(ego), file = paste0("results/GO_enrichment_results_", export_module, ".csv"), row.names = FALSE)
  
  cat("GO analizi bitti!\n")
} else {
  cat("bu modul icin anlamli GO terimi bulunamadi\n")
}
