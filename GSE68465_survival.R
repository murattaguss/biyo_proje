# survival analizi - genlerin hayatta kalmayla iliskisine bakiyoruz
library(survival)
library(survminer)

dir.create("plots", showWarnings = FALSE)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")
dge_sorted <- read.csv("results/GSE68465-dge.csv", row.names = 1)

# hem up hem down genlerden secmemiz lazim, 3 down 2 up aldik
up_genes <- rownames(dge_sorted[dge_sorted$logFC > 0, ])
down_genes <- rownames(dge_sorted[dge_sorted$logFC < 0, ])
top_5_genes <- c(down_genes[1:3], up_genes[1:2])
cat("secilen genler:", paste(top_5_genes, collapse = ", "), "\n")

# her gen icin ayri KM egrisi ciziyoruz
for (gene in top_5_genes) {
  cat("KM egrisi:", gene, "\n")
  
  df_surv <- data.frame(
    expression = as.numeric(GSE68465[, gene]),
    time = df_label$time,
    status = df_label$event
  )
  rownames(df_surv) <- rownames(df_label)
  df_surv <- df_surv[complete.cases(df_surv), ]
  
  # medyana gore high/low diye boluyoruz
  med <- median(df_surv$expression, na.rm = TRUE)
  df_surv$group <- ifelse(df_surv$expression > med, "High", "Low")
  df_surv$group <- factor(df_surv$group, levels = c("Low", "High"))
  
  surv_obj <- Surv(time = df_surv$time, event = df_surv$status)
  fit_km <- survfit(surv_obj ~ group, data = df_surv)
  
  # pval=TRUE log-rank p degerini yazdiriyo, risk.table altta hasta sayisini gosteriyo
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

# en anlamli gen icin cox regresyonu
# sadece KM yetmez, genin etkisi yas/evre/cinsiyet gibi seylerden bagimsiz mi bunu gormeliyiz
top_overall <- rownames(dge_sorted)[1]
cat("cox modeli kuruluyor:", top_overall, "\n")

df_surv_top <- data.frame(
  expression = as.numeric(GSE68465[, top_overall]),
  time = df_label$time,
  status = df_label$event
)
rownames(df_surv_top) <- rownames(df_label)
df_surv_top <- df_surv_top[complete.cases(df_surv_top), ]

med_top <- median(df_surv_top$expression, na.rm = TRUE)
df_surv_top$group <- ifelse(df_surv_top$expression > med_top, "High", "Low")
df_surv_top$group <- factor(df_surv_top$group, levels = c("Low", "High"))

# klinik degiskenleri ekliyoruz
df_surv_top$age <- df_label$age[match(rownames(df_surv_top), rownames(df_label))]
df_surv_top$sex <- df_label$sex[match(rownames(df_surv_top), rownames(df_label))]
df_surv_top$stage <- df_label$binary[match(rownames(df_surv_top), rownames(df_label))]
df_surv_top$smoking <- df_label$smoking[match(rownames(df_surv_top), rownames(df_label))]

cox <- coxph(
  Surv(time, status) ~ group + age + sex + stage + smoking,
  data = df_surv_top
)
print(summary(cox))

# forest plot - hazard ratio degerlerini gosteriyo
cat("forest plot...\n")
png(paste0("plots/forest_", top_overall, ".png"), width=12, height=6, units="in", res=300)
p_forest <- ggforest(cox, data = df_surv_top)
print(p_forest)
dev.off()

cat("survival analizi bitti!\n")
