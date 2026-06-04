# survival analizi
library(survival)
library(survminer)

dir.create("plots", showWarnings = FALSE)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")
dge_sorted <- read.csv("results/GSE68465-dge.csv", row.names = 1)

# up ve down genler
up_genes <- rownames(dge_sorted[dge_sorted$logFC > 0, ])
down_genes <- rownames(dge_sorted[dge_sorted$logFC < 0, ])
top_5_genes <- c(down_genes[1:3], up_genes[1:2])
cat("secilen genler:", paste(top_5_genes, collapse = ", "), "\n")

# KM egrisi
for (gene in top_5_genes) {
  cat("KM egrisi:", gene, "\n")
  
  matched_rows <- match(rownames(GSE68465), rownames(df_label))
  df_surv <- data.frame(
    expression = as.numeric(GSE68465[, gene]),
    time = df_label$time[matched_rows],
    status = df_label$event[matched_rows]
  )
  rownames(df_surv) <- rownames(GSE68465)
  df_surv <- df_surv[complete.cases(df_surv), ]
  
  # mediana gore boluyoruz
  med <- median(df_surv$expression, na.rm = TRUE)
  df_surv$group <- ifelse(df_surv$expression > med, "High", "Low")
  df_surv$group <- factor(df_surv$group, levels = c("Low", "High"))
  
  surv_obj <- Surv(time = df_surv$time, event = df_surv$status)
  fit_km <- survfit(surv_obj ~ group, data = df_surv)
  
  # pval ve risk tablosu
  png(paste0("plots/survival_", gene, ".png"), width=12, height=6, units="in", res=300)
  p_surv <- ggsurvplot(
    fit_km,
    data = df_surv,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    palette = c("#00A087", "#E64B35"),
    title = paste("Kaplan-Meier Survival:", gene),
    legend.labs = c("Low", "High")
  )
  print(p_surv)
  dev.off()
}

# cox regresyonu - top 5 gen + klinik degiskenlerle cok degiskenli model
cat("cox modeli kuruluyor (multivariate)...\n")

matched_cox <- match(rownames(GSE68465), rownames(df_label))
df_cox <- data.frame(
  time = df_label$time[matched_cox],
  status = df_label$event[matched_cox],
  age = df_label$age[matched_cox],
  sex = df_label$sex[matched_cox],
  smoking = df_label$smoking[matched_cox],
  grade = df_label$grade[matched_cox]
)
rownames(df_cox) <- rownames(GSE68465)

# top 5 genin ifadesini medyana gore High/Low yapip ekliyoruz
for (i in seq_along(top_5_genes)) {
  g <- top_5_genes[i]
  vals <- as.numeric(GSE68465[, g])
  med_g <- median(vals, na.rm = TRUE)
  df_cox[[paste0("gene", i)]] <- factor(ifelse(vals > med_g, "High", "Low"), levels = c("Low", "High"))
}

df_cox <- df_cox[complete.cases(df_cox), ]

cox <- coxph(
  Surv(time, status) ~ gene1 + gene2 + gene3 + gene4 + gene5 + age + sex + smoking + grade,
  data = df_cox
)
print(summary(cox))

# forest plot
cat("forest plot...\n")
png("plots/forest_multivariate.png", width=12, height=8, units="in", res=300)
p_forest <- ggforest(cox, data = df_cox)
print(p_forest)
dev.off()

cat("survival analizi bitti!\n")
