# makine ogrenmesi - alive vs dead siniflandirmasi
# 3 farkli ozellik secimi (setup) ve 3 farkli model deniyoruz
library(randomForest)
library(e1071)
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("veriler yukleniyor...\n")
GSE68465 <- readRDS("data/GSE68465.rds")
df_label <- readRDS("data/GSE68465_label.rds")

if (!all(rownames(GSE68465) == rownames(df_label))) {
  df_label <- df_label[match(rownames(GSE68465), rownames(df_label)), , drop = FALSE]
}

vital_binary <- factor(ifelse(df_label$event == 1, "Dead", "Alive"), levels = c("Alive", "Dead"))
n_samples <- nrow(GSE68465)

# --- SETUP 1: biyolojik ozellikler ---
# DGE'den en anlamli genleri + WGCNA hub genlerini aliyoruz
dge_sorted <- read.csv("results/GSE68465-dge.csv", row.names = 1)
up_genes <- rownames(dge_sorted[dge_sorted$logFC > 0, ])
down_genes <- rownames(dge_sorted[dge_sorted$logFC < 0, ])
top_degs <- c(down_genes[1:30], up_genes[1:20])

hubs_df <- read.csv("results/hubsInEachModule.csv", stringsAsFactors = FALSE)
top_hubs <- c()
for (mod in unique(hubs_df$Module)) {
  mod_subset <- hubs_df[hubs_df$Module == mod, ]
  top_hubs <- c(top_hubs, mod_subset$Gene[1])
}

bio_genes <- unique(c(top_degs, top_hubs))
bio_genes <- bio_genes[bio_genes %in% colnames(GSE68465)]
cat("setup-1 gen sayisi:", length(bio_genes), "\n")

# --- SETUP 2 ve 3 icin en degisken 4000 gen ---
gene_vars <- apply(GSE68465, 2, var)
top_4000_genes <- names(sort(gene_vars, decreasing = TRUE)[1:4000])
df_4000 <- GSE68465[, top_4000_genes]

# stratified 5-fold CV
# seed sabitliyoruz ki her seferinde ayni sonuc gelsin
set.seed(42)
folds <- list()
alive_indices <- which(vital_binary == "Alive")
dead_indices <- which(vital_binary == "Dead")

alive_indices <- sample(alive_indices)
dead_indices <- sample(dead_indices)

alive_folds <- split(alive_indices, cut(seq_along(alive_indices), 5, labels = FALSE))
dead_folds <- split(dead_indices, cut(seq_along(dead_indices), 5, labels = FALSE))

for (i in 1:5) {
  folds[[i]] <- c(alive_folds[[i]], dead_folds[[i]])
}

# sonuclari kaydetmek icin bos listeler
models <- c("RF", "SVM", "GLM")
setups <- c("Setup1", "Setup2", "Setup3")

pred_classes <- list()
pred_probs <- list()

for (m in models) {
  for (s in setups) {
    key <- paste0(m, "_", s)
    pred_classes[[key]] <- factor(rep(NA, n_samples), levels = c("Alive", "Dead"))
    pred_probs[[key]] <- rep(NA, n_samples)
  }
}

# --- CV DONGUSU ---
for (fold_idx in 1:5) {
  cat("fold", fold_idx, "...\n")
  test_idx <- folds[[fold_idx]]
  train_idx <- setdiff(1:n_samples, test_idx)
  
  y_train <- vital_binary[train_idx]
  y_test <- vital_binary[test_idx]
  
  # setup 1
  X_train_s1 <- as.matrix(GSE68465[train_idx, bio_genes])
  X_test_s1 <- as.matrix(GSE68465[test_idx, bio_genes])
  
  # setup 2 - LASSO ile ozellik secimi
  # onemli: data leakage olmasin diye sadece train uzerinde yapiyoruz
  X_train_4000 <- as.matrix(df_4000[train_idx, ])
  y_train_num <- ifelse(y_train == "Dead", 1, 0)
  
  cv_lasso <- cv.glmnet(X_train_4000, y_train_num, family = "binomial", alpha = 1)
  coefs <- coef(cv_lasso, s = "lambda.min")
  selected_genes <- rownames(coefs)[which(coefs[, 1] != 0)]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]
  
  # hic gen secilmezse en korele 10 geni al
  if (length(selected_genes) == 0) {
    cors <- abs(cor(X_train_4000, y_train_num))
    selected_genes <- rownames(cors)[order(-cors)[1:10]]
  }
  # cok fazla gen secilirse 30'a kisitla, yoksa glm cokuyo
  if (length(selected_genes) > 30) {
    lasso_coefs <- coefs[selected_genes, 1]
    selected_genes <- names(sort(abs(lasso_coefs), decreasing = TRUE)[1:30])
  }
  
  X_train_s2 <- as.matrix(df_4000[train_idx, selected_genes, drop = FALSE])
  X_test_s2 <- as.matrix(df_4000[test_idx, selected_genes, drop = FALSE])
  
  # setup 3 - PCA ile boyut indirgeme
  # bunu da sadece train uzerinde egitiyoruz
  pca_train <- prcomp(df_4000[train_idx, ], scale. = TRUE)
  X_train_s3 <- pca_train$x[, 1:10]
  X_test_s3 <- predict(pca_train, newdata = df_4000[test_idx, ])[, 1:10]
  
  for (s in setups) {
    if (s == "Setup1") {
      X_tr <- X_train_s1; X_te <- X_test_s1
    } else if (s == "Setup2") {
      X_tr <- X_train_s2; X_te <- X_test_s2
    } else {
      X_tr <- X_train_s3; X_te <- X_test_s3
    }
    
    # random forest
    fit_rf <- randomForest(x = X_tr, y = y_train, ntree = 500)
    pred_classes[[paste0("RF_", s)]][test_idx] <- predict(fit_rf, X_te, type = "response")
    pred_probs[[paste0("RF_", s)]][test_idx] <- predict(fit_rf, X_te, type = "prob")[, "Dead"]
    
    # svm
    fit_svm <- svm(x = X_tr, y = y_train, probability = TRUE)
    svm_pred <- predict(fit_svm, X_te, probability = TRUE)
    pred_classes[[paste0("SVM_", s)]][test_idx] <- svm_pred
    pred_probs[[paste0("SVM_", s)]][test_idx] <- attr(svm_pred, "probabilities")[, "Dead"]
    
    # lojistik regresyon
    df_tr <- as.data.frame(X_tr)
    df_tr$label <- y_train
    fit_glm <- glm(label ~ ., data = df_tr, family = binomial)
    
    df_te <- as.data.frame(X_te)
    probs_glm <- predict(fit_glm, newdata = df_te, type = "response")
    pred_probs[[paste0("GLM_", s)]][test_idx] <- probs_glm
    pred_classes[[paste0("GLM_", s)]][test_idx] <- factor(ifelse(probs_glm > 0.5, "Dead", "Alive"), levels = c("Alive", "Dead"))
  }
}

# performans metrikleri
calc_metrics <- function(true_y, pred_y, pred_p) {
  conf <- table(True = true_y, Pred = pred_y)
  TN <- conf["Alive", "Alive"]
  FN <- conf["Dead", "Alive"]
  FP <- conf["Alive", "Dead"]
  TP <- conf["Dead", "Dead"]
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  roc_obj <- roc(true_y, pred_p, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  if (is.nan(precision)) precision <- 0
  if (is.nan(f1)) f1 <- 0
  
  return(c(Accuracy = accuracy, Precision = precision, Recall = recall, F1 = f1, AUC = auc_val))
}

results_list <- list()
for (m in models) {
  for (s in setups) {
    key <- paste0(m, "_", s)
    metrics <- calc_metrics(vital_binary, pred_classes[[key]], pred_probs[[key]])
    results_list[[key]] <- data.frame(
      Model = m, Setup = s,
      Accuracy = metrics["Accuracy"], Precision = metrics["Precision"],
      Recall = metrics["Recall"], F1_score = metrics["F1"], AUC = metrics["AUC"],
      stringsAsFactors = FALSE
    )
  }
}

df_perf <- do.call(rbind, results_list)
rownames(df_perf) <- NULL
write.csv(df_perf, "results/ML_performance_comparison.csv", row.names = FALSE)
print(df_perf)

# ROC egrileri - her setup icin ayri grafik
cat("ROC grafikleri ciziliyor...\n")
png("plots/ROC_comparison.png", width=15, height=5, units="in", res=300)
par(mfrow = c(1, 3))

colors <- c("RF" = "#377EB8", "SVM" = "#E41A1C", "GLM" = "#4DAF4A")

for (s in setups) {
  plot(0, 0, type="n", xlim=c(1, 0), ylim=c(0, 1),
       xlab="Specificity", ylab="Sensitivity",
       main=paste("ROC Curves -", s))
  abline(a=1, b=-1, lty=2, col="grey")
  
  for (m in models) {
    key <- paste0(m, "_", s)
    roc_obj <- roc(vital_binary, pred_probs[[key]], quiet = TRUE)
    plot(roc_obj, add=TRUE, col=colors[m], lwd=2)
  }
  
  legend_labels <- sapply(models, function(m) {
    key <- paste0(m, "_", s)
    roc_obj <- roc(vital_binary, pred_probs[[key]], quiet = TRUE)
    paste0(m, " (AUC=", round(as.numeric(auc(roc_obj)), 3), ")")
  })
  legend("bottomright", legend=legend_labels, col=colors, lwd=2, cex=0.9)
}
dev.off()

# en iyi modelin confusion matrix'i
best_idx <- which.max(df_perf$F1_score)
best_model <- df_perf$Model[best_idx]
best_setup <- df_perf$Setup[best_idx]
best_key <- paste0(best_model, "_", best_setup)

cat("en iyi model:", best_model, "-", best_setup, "\n")

conf_mat <- table(True = vital_binary, Pred = pred_classes[[best_key]])
conf_df <- as.data.frame(conf_mat)

png("plots/confusion_matrix_best.png", width=6, height=5, units="in", res=300)
p_cm <- ggplot(conf_df, aes(x = Pred, y = True, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 6, color = "black") +
  scale_fill_gradient(low = "#E0F2F1", high = "#00796B") +
  theme_minimal() +
  labs(
    title = paste("Confusion Matrix -", best_model, "(", best_setup, ")"),
    x = "Predicted", y = "True"
  ) +
  theme(axis.text = element_text(size = 12),
        title = element_text(size = 10, face = "bold"))
print(p_cm)
dev.off()

cat("ML bitti!\n")
