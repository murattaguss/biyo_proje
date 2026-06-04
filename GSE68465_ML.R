# makine ogrenmesi
# 3 setup ve 3 model deniyoruz
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

grade_binary <- factor(df_label$grade, levels = c("Well", "Poorly"))
keep_samples <- !is.na(grade_binary)
GSE68465 <- GSE68465[keep_samples, , drop = FALSE]
df_label <- df_label[keep_samples, , drop = FALSE]
grade_binary <- droplevels(factor(df_label$grade, levels = c("Well", "Poorly")))
n_samples <- nrow(GSE68465)

# setup 1
# DGE ve WGCNA genleri
dge_sorted <- read.csv("results/GSE68465-dge.csv", row.names = 1)
up_genes <- rownames(dge_sorted[dge_sorted$logFC > 0, ])
down_genes <- rownames(dge_sorted[dge_sorted$logFC < 0, ])
top_degs <- na.omit(c(head(down_genes, 30), head(up_genes, 20)))

hubs_df <- read.csv("results/hubsInEachModule.csv", stringsAsFactors = FALSE)
top_hubs <- c()
for (mod in unique(hubs_df$Module)) {
  mod_subset <- hubs_df[hubs_df$Module == mod, ]
  top_hubs <- c(top_hubs, mod_subset$Gene[1])
}

bio_genes <- unique(c(top_degs, top_hubs))
bio_genes <- bio_genes[bio_genes %in% colnames(GSE68465)]
cat("setup-1 gen sayisi:", length(bio_genes), "\n")

# setup 2 ve 3
df_all <- GSE68465

# 5-fold CV
# seed sabit
set.seed(42)
folds <- list()
well_indices <- which(grade_binary == "Well")
poorly_indices <- which(grade_binary == "Poorly")

well_indices <- sample(well_indices)
poorly_indices <- sample(poorly_indices)

well_folds <- split(well_indices, cut(seq_along(well_indices), 5, labels = FALSE))
poorly_folds <- split(poorly_indices, cut(seq_along(poorly_indices), 5, labels = FALSE))

for (i in 1:5) {
  folds[[i]] <- c(well_folds[[i]], poorly_folds[[i]])
}

# sonuc listeleri
models <- c("RF", "SVM", "GLM")
setups <- c("Setup1", "Setup2", "Setup3")

pred_classes <- list()
pred_probs <- list()

for (m in models) {
  for (s in setups) {
    key <- paste0(m, "_", s)
    pred_classes[[key]] <- factor(rep(NA, n_samples), levels = c("Well", "Poorly"))
    pred_probs[[key]] <- rep(NA, n_samples)
  }
}

# CV dongusu
for (fold_idx in 1:5) {
  cat("fold", fold_idx, "...\n")
  test_idx <- folds[[fold_idx]]
  train_idx <- setdiff(1:n_samples, test_idx)
  
  y_train <- grade_binary[train_idx]
  y_test <- grade_binary[test_idx]
  
  # setup 1
  X_train_s1 <- as.matrix(GSE68465[train_idx, bio_genes])
  X_test_s1 <- as.matrix(GSE68465[test_idx, bio_genes])
  
  # setup 2 - LASSO
  # data leakage olmasin diye train'de calisiyoruz
  X_train_all <- as.matrix(df_all[train_idx, ])
  y_train_num <- ifelse(y_train == "Poorly", 1, 0)
  
  cv_lasso <- cv.glmnet(X_train_all, y_train_num, family = "binomial", alpha = 1)
  coefs <- coef(cv_lasso, s = "lambda.min")
  selected_genes <- rownames(coefs)[which(coefs[, 1] != 0)]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]
  
  # gen secilmezse en korelileri al
  if (length(selected_genes) == 0) {
    cors <- abs(cor(X_train_all, y_train_num))
    selected_genes <- rownames(cors)[order(-cors)[1:10]]
  }
  # cok gen secilirse kisitla
  if (length(selected_genes) > 30) {
    lasso_coefs <- coefs[selected_genes, 1]
    selected_genes <- names(sort(abs(lasso_coefs), decreasing = TRUE)[1:30])
  }
  
  X_train_s2 <- as.matrix(df_all[train_idx, selected_genes, drop = FALSE])
  X_test_s2 <- as.matrix(df_all[test_idx, selected_genes, drop = FALSE])
  
  # setup 3 - PCA
  # train'de egitiyoruz
  pca_train <- prcomp(df_all[train_idx, ], scale. = TRUE)
  X_train_s3 <- pca_train$x[, 1:10]
  X_test_s3 <- predict(pca_train, newdata = df_all[test_idx, ])[, 1:10]
  
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
    pred_probs[[paste0("RF_", s)]][test_idx] <- predict(fit_rf, X_te, type = "prob")[, "Poorly"]
    
    # svm
    fit_svm <- svm(x = X_tr, y = y_train, probability = TRUE)
    svm_pred <- predict(fit_svm, X_te, probability = TRUE)
    pred_classes[[paste0("SVM_", s)]][test_idx] <- svm_pred
    probs_svm <- attr(svm_pred, "probabilities")
    pred_probs[[paste0("SVM_", s)]][test_idx] <- probs_svm[, which(colnames(probs_svm) == "Poorly")]
    
    # lojistik regresyon
    df_tr <- as.data.frame(X_tr)
    df_tr$label <- y_train
    fit_glm <- glm(label ~ ., data = df_tr, family = binomial)
    
    df_te <- as.data.frame(X_te)
    probs_glm <- predict(fit_glm, newdata = df_te, type = "response")
    pred_probs[[paste0("GLM_", s)]][test_idx] <- probs_glm
    pred_classes[[paste0("GLM_", s)]][test_idx] <- factor(ifelse(probs_glm > 0.5, "Poorly", "Well"), levels = c("Well", "Poorly"))
  }
}

# metrikler
calc_metrics <- function(true_y, pred_y, pred_p) {
  conf <- table(factor(true_y, levels = c("Well", "Poorly")),
                factor(pred_y, levels = c("Well", "Poorly")))
  TN <- conf["Well", "Well"]
  FN <- conf["Poorly", "Well"]
  FP <- conf["Well", "Poorly"]
  TP <- conf["Poorly", "Poorly"]
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  roc_obj <- roc(response = true_y, predictor = pred_p, levels = c("Well", "Poorly"), direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  if (is.nan(precision)) precision <- 0
  if (is.nan(f1)) f1 <- 0
  
  return(c(Accuracy = accuracy, Precision = precision, Recall = recall, F1 = f1, AUC = auc_val))
}

results_list <- list()
for (m in models) {
  for (s in setups) {
    key <- paste0(m, "_", s)
    metrics <- calc_metrics(grade_binary, pred_classes[[key]], pred_probs[[key]])
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

# ROC egrileri
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
    roc_obj <- roc(response = grade_binary, predictor = pred_probs[[key]], levels = c("Well", "Poorly"), direction = "<", quiet = TRUE)
    plot(roc_obj, add=TRUE, col=colors[m], lwd=2)
  }
  
  legend_labels <- sapply(models, function(m) {
    key <- paste0(m, "_", s)
    roc_obj <- roc(response = grade_binary, predictor = pred_probs[[key]], levels = c("Well", "Poorly"), direction = "<", quiet = TRUE)
    paste0(m, " (AUC=", round(as.numeric(auc(roc_obj)), 3), ")")
  })
  legend("bottomright", legend=legend_labels, col=colors, lwd=2, cex=0.9)
}
dev.off()

# tum modeller icin confusion matrix (3x3 grid)
cm_plots <- list()
for (m in models) {
  for (s in setups) {
    key <- paste0(m, "_", s)
    conf_mat <- table(True = grade_binary, Pred = pred_classes[[key]])
    conf_df <- as.data.frame(conf_mat)
    
    p_cm <- ggplot(conf_df, aes(x = Pred, y = True, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), size = 5, color = "black") +
      scale_fill_gradient(low = "#E0F2F1", high = "#00796B") +
      theme_minimal() +
      labs(title = paste(m, "-", s), x = "Predicted", y = "True") +
      theme(axis.text = element_text(size = 10),
            title = element_text(size = 9, face = "bold"))
    cm_plots[[key]] <- p_cm
  }
}

png("plots/confusion_matrices_all.png", width=15, height=12, units="in", res=300)
gridExtra::grid.arrange(grobs = cm_plots, ncol = 3)
dev.off()

cat("ML bitti!\n")
