# grafikler
library(ggplot2)
library(reshape2)
library(limma)
library(Biobase)
library(ggpubr)

dir.create("plots", showWarnings = FALSE)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")
dge_sorted <- read.csv("results/GSE68465-dge.csv", row.names = 1)

# grade etiketi
df_label$grade_binary <- factor(df_label$grade, levels = c("Well", "Poorly"))
keep_samples <- !is.na(df_label$grade_binary)
GSE68465 <- GSE68465[keep_samples, , drop = FALSE]
df_label <- df_label[keep_samples, , drop = FALSE]
df_label$grade_binary <- droplevels(factor(df_label$grade, levels = c("Well", "Poorly")))

# volcano icin transpose
df_mat <- as.matrix(GSE68465)
mode(df_mat) <- "numeric"
df_mat <- t(df_mat)

eset <- ExpressionSet(assayData = df_mat, phenoData = AnnotatedDataFrame(df_label))
design <- model.matrix(~0 + grade_binary, data = pData(eset))
colnames(design) <- c("grade_binaryWell", "grade_binaryPoorly")
cm <- makeContrasts(PoorlyvWell = grade_binaryPoorly - grade_binaryWell, levels = design)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, number = Inf)

# renk grubu
tt$threshold <- "Not Sig"
tt$threshold[tt$adj.P.Val < 0.05 & tt$logFC > 1]  <- "Up"
tt$threshold[tt$adj.P.Val < 0.05 & tt$logFC < -1] <- "Down"
tt$threshold <- factor(tt$threshold, levels = c("Down", "Not Sig", "Up"))

# volcano
cat("volcano plot...\n")
png("plots/volcanoGSE68465.png", width=12, height=6, units="in", res=300)
p_volcano <- ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot (Poorly vs Well)",
       x = "log2 Fold Change",
       y = "-log10(FDR)")
print(p_volcano)
dev.off()

top_gene <- rownames(dge_sorted)[1]
cat("en anlamli gen:", top_gene, "\n")

df_plot <- data.frame(
  expression = as.numeric(GSE68465[, top_gene]),
  group = df_label$grade_binary
)
df_plot$group <- factor(df_plot$group, levels = c("Well", "Poorly"))

# boxplot
cat("boxplot...\n")
png(paste0("plots/boxplot_", top_gene, ".png"), width=12, height=6, units="in", res=300)
p_box <- ggplot(df_plot, aes(x = group, y = expression, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5) +
  scale_fill_manual(values = c("#00A087", "#E64B35")) +
  theme_minimal() +
  labs(title = paste("Expression of", top_gene, "by Histologic Grade"),
       x = "Group",
       y = "Expression")
print(p_box)
dev.off()

# violin
cat("violin...\n")
png(paste0("plots/violin_", top_gene, ".png"), width=12, height=6, units="in", res=300)
p_violin <- ggplot(df_plot, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5) +
  scale_fill_manual(values = c("#00A087", "#E64B35")) +
  theme_minimal() +
  labs(title = paste("Expression of", top_gene, "by Histologic Grade"),
       x = "Group",
       y = "Expression")
print(p_violin)
dev.off()

# top 5
top_genes <- rownames(dge_sorted)[1:5]
df_top <- GSE68465[, top_genes, drop = FALSE]
df_top$sample <- rownames(df_top)

# uzun form
df_long <- melt(
  df_top,
  id.vars = "sample",
  variable.name = "gene",
  value.name = "expression"
)
df_long$group <- df_label[df_long$sample, "grade_binary"]
df_long$expression <- as.numeric(df_long$expression)
df_long$group <- factor(df_long$group, levels = c("Well", "Poorly"))

# top 5 violin
cat("top 5 violin...\n")
png("plots/violinTop5.png", width=12, height=6, units="in", res=300)
p_top5_violin <- ggplot(df_long, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  facet_wrap(~gene, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 1.5) +
  scale_fill_manual(values = c("#00A087", "#E64B35")) +
  theme_minimal() +
  labs(title = "Top 5 DEGs - Well vs Poorly", x = "Group", y = "Expression")
print(p_top5_violin)
dev.off()

# top 5 boxplot
cat("top 5 boxplot...\n")
png("plots/boxplotTop5.png", width=12, height=6, units="in", res=300)
p_top5_box <- ggplot(df_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  facet_wrap(~gene, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label.x = 1.5) +
  scale_fill_manual(values = c("#00A087", "#E64B35")) +
  theme_minimal() +
  labs(title = "Top 5 DEGs Boxplot - Well vs Poorly", x = "Group", y = "Expression")
print(p_top5_box)
dev.off()

cat("grafikler bitti!\n")
