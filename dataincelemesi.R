library(GEOquery)

gse <- getGEO("GSE68465", GSEMatrix = TRUE)
df_label_raw <- pData(gse[[1]])
print(colnames(df_label_raw))
print(table(df_label_raw$`disease_state:ch1`, useNA = "always"))
print(table(df_label_raw$`histologic_grade:ch1`, useNA = "always"))
print(table(df_label_raw$`vital_status:ch1`, useNA = "always"))
print(head(df_label_raw$`histologic_grade:ch1`, 20))
