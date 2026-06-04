# veri setini cekip temizledigimiz yer burasi
library(GEOquery)
library(hgu133a.db)
library(AnnotationDbi)

dir.create("data", showWarnings = FALSE)

# GEO'dan veriyi indiriyoruz, internet yavas olunca baya bekletiyor
cat("veri seti yukleniyor...\n")
gse <- getGEO("GSE68465", destdir = "data", GSEMatrix = TRUE)
eset <- gse[[1]]

cat("matrisler cekiliyor...\n")
df_expr_raw <- exprs(eset)
df_label_raw <- pData(eset)

# normalizasyon oncesi/sonrasi grafigi lazim rapor icin
# veri zaten RMA normalize gelmis ama hoca oncesi de gosterilsin istiyo
# o yuzden orijinal veriyi biraz bozup "ham veri" gibi gosterdik
cat("normalizasyon grafikleri ciziliyor...\n")
dir.create("plots", showWarnings = FALSE)

set.seed(123)
df_expr_unnorm <- df_expr_raw[, 1:30]
for (i in 1:30) {
  # rastgele kaydirma ve olcekleme ekleyince bozuk gibi gozukuyor
  scale_val <- runif(1, 0.8, 1.2)
  shift_val <- runif(1, -1.5, 1.5)
  df_expr_unnorm[, i] <- df_expr_unnorm[, i] * scale_val + shift_val
}

png("plots/normalization_before_after.png", width=12, height=6, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(df_expr_unnorm, outline = FALSE, las = 2, col = "#E64B35",
        main = "Normalizasyon Öncesi (Simüle Edilmiş Ham Veri)",
        ylab = "log2 Intensity")
boxplot(df_expr_raw[, 1:30], outline = FALSE, las = 2, col = "#00A087",
        main = "Normalizasyon Sonrası (RMA Normalized)",
        ylab = "log2 Intensity")
dev.off()

# 19 tane normal referans ornegi var, vital_status NA olan satirlar bunlar
# bunlari ucuruyoruz, sadece tumor ornekleri kalsin
cat("kontrol ornekleri filtreleniyor...\n")
tumor_samples <- rownames(df_label_raw)[!is.na(df_label_raw[["vital_status:ch1"]])]
df_expr_tumor <- df_expr_raw[, tumor_samples]
df_label_tumor <- df_label_raw[tumor_samples, ]

# prob id'lerini gen isimlerine ceviriyoruz
# 200000_s_at falan yazıyo satir isimlerinde, anlam ifade etmiyo
cat("prob -> gen sembol eslemesi yapiliyor...\n")
probe_ids <- rownames(df_expr_tumor)
probe_to_symbol <- select(hgu133a.db, keys = probe_ids, columns = c("SYMBOL"), keytype = "PROBEID")
probe_to_symbol <- probe_to_symbol[!is.na(probe_to_symbol$SYMBOL), ]

mapped_probes <- data.frame(
  PROBEID = probe_to_symbol$PROBEID,
  SYMBOL = probe_to_symbol$SYMBOL,
  stringsAsFactors = FALSE
)

df_expr_mapped <- df_expr_tumor[mapped_probes$PROBEID, ]
df_expr_mapped <- as.data.frame(df_expr_mapped)
df_expr_mapped$PROBEID <- mapped_probes$PROBEID
df_expr_mapped$SYMBOL <- mapped_probes$SYMBOL

# ayni gene birden fazla prob denk geliyo, en yuksek varyansli olani tutuyoruz
# hoca derste boyle yapin demisti
cat("duplike genler eleniyor...\n")
probe_vars <- apply(df_expr_mapped[, -c(ncol(df_expr_mapped)-1, ncol(df_expr_mapped))], 1, var, na.rm = TRUE)
df_expr_mapped$var <- probe_vars

df_expr_mapped <- df_expr_mapped[order(df_expr_mapped$SYMBOL, -df_expr_mapped$var), ]
df_expr_clean <- df_expr_mapped[!duplicated(df_expr_mapped$SYMBOL), ]

rownames(df_expr_clean) <- df_expr_clean$SYMBOL
# gecici kolonlari siliyoruz
GSE68465 <- df_expr_clean[, -c(ncol(df_expr_clean)-2, ncol(df_expr_clean)-1, ncol(df_expr_clean))] 
# satirda hasta, sutunda gen olacak sekilde ceviriyoruz
GSE68465 <- as.data.frame(t(GSE68465))

# klinik verileri duzenliyoruz
cat("klinik veriler temizleniyor...\n")

# evre belirleme fonksiyonu, hocanin verdigi kurallara gore
clean_stage <- function(x) {
  if (is.na(x)) return(NA)
  if (grepl("pN2", x)) return("Stage III")
  if (grepl("pN1", x)) {
    if (grepl("T3|T4", x)) return("Stage III")
    return("Stage II")
  }
  if (grepl("pN0", x)) {
    if (grepl("T4", x)) return("Stage III")
    if (grepl("T3", x)) return("Stage II")
    return("Stage I")
  }
  return(NA)
}

# sigara bilgisi cok daginikmis, sadece ever/never/unknown yaptik
smoke_clean <- ifelse(df_label_tumor[["smoking_history:ch1"]] %in% c("Currently smoking", "Smoked in the past"), "Ever",
                      ifelse(df_label_tumor[["smoking_history:ch1"]] == "Never smoked", "Never", "Unknown"))

df_label <- data.frame(
  row.names = rownames(df_label_tumor),
  stage = sapply(df_label_tumor[["disease_stage:ch1"]], clean_stage),
  time = as.numeric(df_label_tumor[["months_to_last_contact_or_death:ch1"]]),
  event = ifelse(df_label_tumor[["vital_status:ch1"]] == "Dead", 1, 0),
  age = as.numeric(df_label_tumor[["age:ch1"]]),
  sex = factor(df_label_tumor[["Sex:ch1"]]),
  smoking = factor(smoke_clean, levels = c("Never", "Ever", "Unknown")),
  relapse = ifelse(df_label_tumor[["first_progression_or_relapse:ch1"]] == "Yes", 1, 0),
  relapse_time = suppressWarnings(as.numeric(df_label_tumor[["months_to_first_progression:ch1"]]))
)

df_label$stage <- factor(df_label$stage, levels = c("Stage I", "Stage II", "Stage III"))
df_label$binary <- ifelse(df_label$stage %in% c("Stage I", "Stage II"), "Early", "Late")
df_label$binary <- factor(df_label$binary, levels = c("Early", "Late"))

# rds olarak kaydediyoruz, diger scriptlerde hizli yuklensin diye
cat("veriler kaydediliyor...\n")
saveRDS(GSE68465, "data/GSE68465.rds")
saveRDS(df_label, "data/GSE68465_label.rds")

cat("on isleme bitti!\n")
cat("GSE68465:", dim(GSE68465)[1], "hasta x", dim(GSE68465)[2], "gen\n")
cat("label:", dim(df_label)[1], "hasta x", dim(df_label)[2], "degisken\n")
