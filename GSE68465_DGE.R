# DGE analizi - olenler vs yasayanlar arasinda hangi genler farkli ifade ediliyo
library(limma)
library(tidyr)
library(Biobase)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")

# siralamalarin ayni oldugundan emin olalim yoksa her sey cop olur
cat("siralama kontrolu...\n")
if (!all(rownames(GSE68465) == rownames(df_label))) {
  cat("siralama kayik, duzeltiliyor...\n")
  df_label <- df_label[match(rownames(GSE68465), rownames(df_label)), , drop = FALSE]
}
print(all(rownames(GSE68465) == rownames(df_label)))

# limma genleri satirda istiyo, biz sutuna koymusduk, o yuzden t() aliyoruz
df_mat <- as.matrix(GSE68465)
mode(df_mat) <- "numeric"
df_mat <- t(df_mat) 

# 1 = olmus, 0 = yasıyo. faktor yapalim
df_label$vital_binary <- ifelse(df_label$event == 1, "Dead", "Alive")
df_label$vital_binary <- factor(df_label$vital_binary, levels = c("Alive", "Dead"))

cat("grup dagilimi:\n")
print(table(df_label$vital_binary))

eset <- ExpressionSet(assayData = df_mat, phenoData = AnnotatedDataFrame(df_label))

# ~0 koyunca intercept olmuyo, gruplari direkt kiyasliyoruz
design <- model.matrix(~0 + vital_binary, data = pData(eset))
colnames(design) <- c("vital_binaryAlive", "vital_binaryDead")
head(design)

# dead - alive kontrastı, yani olenlerde nelerin degistigini ariyoruz
cm <- makeContrasts(AlivevDead = vital_binaryDead - vital_binaryAlive, levels = design)

# limma pipeline: fit -> kontrast uygula -> ebayes duzeltmesi
cat("model fit ediliyor...\n")
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)

# kac gen anlamli cikmis bakalim
results <- decideTests(fit2)
cat("\nanlamli gen sayilari (FDR < 0.05):\n")
print(summary(results))

# butun genlerin tablosu
tt <- topTable(fit2, number = Inf)

# sadece anlamli olanlari filtrele ve logFC'ye gore sirala
sig_genes <- rownames(results)[results[, 1] != 0]
dge_table <- tt[rownames(tt) %in% sig_genes, ]
dge_sorted <- dge_table[order(-abs(dge_table$logFC)), ]

dir.create("results", showWarnings = FALSE)

cat("sonuclar kaydediliyor...\n")
write.csv(dge_sorted, "results/GSE68465-dge.csv", row.names = TRUE)

cat("DGE bitti!\n")
