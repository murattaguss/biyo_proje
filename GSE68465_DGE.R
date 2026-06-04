# DGE analizi - Well vs Poorly
library(limma)
library(tidyr)
library(Biobase)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")

# siralamayi kontrol et
cat("siralama kontrolu...\n")
if (!all(rownames(GSE68465) == rownames(df_label))) {
  cat("siralama kayik, duzeltiliyor...\n")
  df_label <- df_label[match(rownames(GSE68465), rownames(df_label)), , drop = FALSE]
}
print(all(rownames(GSE68465) == rownames(df_label)))

# grade etiketi
df_label$grade_binary <- factor(df_label$grade, levels = c("Well", "Poorly"))
keep_samples <- !is.na(df_label$grade_binary)
GSE68465 <- GSE68465[keep_samples, , drop = FALSE]
df_label <- df_label[keep_samples, , drop = FALSE]

# limma icin transpose
df_mat <- as.matrix(GSE68465)
mode(df_mat) <- "numeric"
df_mat <- t(df_mat)

cat("grup dagilimi:\n")
print(table(df_label$grade_binary))

eset <- ExpressionSet(assayData = df_mat, phenoData = AnnotatedDataFrame(df_label))

# intercept olmasin
design <- model.matrix(~0 + grade_binary, data = pData(eset))
colnames(design) <- c("grade_binaryWell", "grade_binaryPoorly")
head(design)

# poorly - well
cm <- makeContrasts(PoorlyvWell = grade_binaryPoorly - grade_binaryWell, levels = design)

# limma
cat("model fit ediliyor...\n")
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)

# anlamli gen sayisi
results <- decideTests(fit2, p.value = 0.05, lfc = 1)
cat("\nanlamli gen sayilari (FDR < 0.05):\n")
print(summary(results))

# tum genler
tt <- topTable(fit2, number = Inf)

# filtrele ve sirala
sig_genes <- rownames(results)[results[, 1] != 0]
dge_table <- tt[rownames(tt) %in% sig_genes, ]
dge_sorted <- dge_table[order(-abs(dge_table$logFC)), ]

dir.create("results", showWarnings = FALSE)

cat("sonuclar kaydediliyor...\n")
write.csv(dge_sorted, "results/GSE68465-dge.csv", row.names = TRUE)

cat("DGE bitti!\n")
