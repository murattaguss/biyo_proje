library(GEOquery)
library(affy)
library(hgu133a.db)
library(AnnotationDbi)

dir.create("data/raw", showWarnings=F, recursive=T)
dir.create("plots", showWarnings=F)

# metadatayi cekelim
gse <- getGEO("GSE68465", destdir="data", GSEMatrix=T)
eset_meta <- gse[[1]]

# raw dosyalari bul ve rma yap
raw_tar <- list.files("data", pattern="GSE68465_RAW\\.tar$", full.names=T, recursive=T)
if(length(raw_tar) == 0) stop("data klasorune GSE68465_RAW.tar at")

untar(raw_tar[1], exdir="data/raw")
cel_files <- list.files("data/raw", pattern="[.]CEL([.]gz)?$", full.names=T, recursive=T, ignore.case=T)
if(length(cel_files) == 0) stop("cel dosyalari bulunamadi")

raw_affy <- ReadAffy(filenames=cel_files)
eset <- rma(raw_affy)

df_expr_raw <- exprs(eset)
df_label_raw <- pData(eset_meta)

# normalizasyon boxplot
n_plot <- min(30, ncol(df_expr_raw))
plot_idx <- seq_len(n_plot)
raw_affy_30 <- raw_affy[, plot_idx]
png("plots/normalization_before_after.png", width=12, height=6, units="in", res=300)
par(mfrow=c(1,2), mar=c(8,4,4,2))
boxplot(log2(intensity(raw_affy_30)), outline=F, las=2, col="#E64B35", main="Norm Oncesi (log2 Raw)", ylab="log2 Intensity")
boxplot(exprs(rma(raw_affy_30)), outline=F, las=2, col="#00A087", main="Norm Sonrasi (RMA)", ylab="log2 Intensity")
dev.off()

# sample isimlerini esle
raw_names <- sampleNames(raw_affy)
gsm_ids <- regmatches(raw_names, regexpr("GSM[0-9]+", raw_names))
sample_map <- data.frame(
  raw_sample = raw_names,
  geo_sample = gsm_ids,
  stringsAsFactors = F
)

matched <- match(sample_map$geo_sample, rownames(df_label_raw))
keep <- !is.na(matched)
sample_map <- sample_map[keep, ]
df_label_raw <- df_label_raw[matched[keep], , drop=F]
rownames(df_label_raw) <- sample_map$raw_sample
df_expr_raw <- df_expr_raw[, sample_map$raw_sample, drop=F]

# normal dokulari (vital_status NA olanlar) atiyoruz sadece tumor kalsin
is_tumor <- !is.na(df_label_raw[["vital_status:ch1"]])
df_expr_tumor <- df_expr_raw[, is_tumor, drop=F]
df_label_tumor <- df_label_raw[is_tumor, , drop=F]

# anotasyon kismi (prob -> gen)
probe_ids <- rownames(df_expr_tumor)
probe_to_sym <- select(hgu133a.db, keys=probe_ids, columns="SYMBOL", keytype="PROBEID")
probe_to_sym <- probe_to_sym[!is.na(probe_to_sym$SYMBOL), ]

df_mapped <- df_expr_tumor[probe_to_sym$PROBEID, ]
df_mapped <- as.data.frame(df_mapped)
df_mapped$SYMBOL <- probe_to_sym$SYMBOL

# ayni gen birden fazla proba denk geliyorsa varyansi en yuksek olani tutalim
df_mapped$var <- apply(df_mapped[, -ncol(df_mapped)], 1, var, na.rm=T)
df_mapped <- df_mapped[order(df_mapped$SYMBOL, -df_mapped$var), ]
df_mapped <- df_mapped[!duplicated(df_mapped$SYMBOL), ]
rownames(df_mapped) <- df_mapped$SYMBOL
df_mapped <- df_mapped[, !(colnames(df_mapped) %in% c("SYMBOL", "var"))]

# matrisi transpose edip df yapiyoruz (hasta satirlarda genler sutunlarda)
GSE68465 <- as.data.frame(t(df_mapped), check.names = FALSE)

# projede istenen top-3000 hvg seciyoruz
gene_var <- apply(GSE68465, 2, var, na.rm=T)
n_hvg <- min(3000, length(gene_var))
top_hvg <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_hvg)]
top_hvg <- intersect(top_hvg, colnames(GSE68465))
GSE68465 <- GSE68465[, top_hvg, drop=F]

# klinik verileri toparlayalim
grade_raw <- toupper(trimws(df_label_tumor[["histologic_grade:ch1"]]))
grade_clean <- ifelse(grade_raw %in% c("--", ""), NA, grade_raw)
grade_clean <- ifelse(
  grade_clean == "WELL DIFFERENTIATED", "Well",
  ifelse(
    grade_clean == "MODERATE DIFFERENTIATION", "Moderate",
    ifelse(grade_clean == "POORLY DIFFERENTIATED", "Poorly", NA)
  )
)

smoke_clean <- ifelse(
  df_label_tumor[["smoking_history:ch1"]] %in% c("Currently smoking", "Smoked in the past"), "Ever",
  ifelse(df_label_tumor[["smoking_history:ch1"]] == "Never smoked", "Never", "Unknown")
)

df_label <- data.frame(
  row.names = rownames(df_label_tumor),
  time = as.numeric(df_label_tumor[["months_to_last_contact_or_death:ch1"]]),
  event = ifelse(df_label_tumor[["vital_status:ch1"]] == "Dead", 1, 0),
  grade = factor(grade_clean, levels = c("Well", "Moderate", "Poorly")),
  age = as.numeric(df_label_tumor[["age:ch1"]]),
  sex = factor(df_label_tumor[["Sex:ch1"]]),
  smoking = factor(smoke_clean, levels=c("Never", "Ever", "Unknown"))
)

# grade verisi eksik olanlari ucur
valid <- !is.na(df_label$grade)
GSE68465 <- GSE68465[valid, , drop=F]
df_label <- df_label[valid, , drop=F]

# diger kodlarda kullanmak icin rds kaydediyoruz
saveRDS(GSE68465, "data/GSE68465.rds")
saveRDS(df_label, "data/GSE68465_label.rds")

cat("veri on isleme bitti!\n")
cat("GSE68465:", nrow(GSE684
