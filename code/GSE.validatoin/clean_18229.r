# ============================================ #
#               clean_18229.R                  #
# ============================================ #

### ----------------- read data -------------------- ###
series <- c("885", "887", "1390", "1708", "5325", "6607", "7504")

### read expression data
for (i in 1:length(series)) {
      ## read expression and annotation data      
      exp.temp <- read.table(paste("data/GSE/GSE18229-GPL", series[i], ".txt", sep = ""), stringsAsFactors = FALSE, header = TRUE, comment.char = "!", row.names = 1, fill = TRUE)
      anno.temp <- read.table(paste("data/GSE/GPL", series[i], ".anno.txt", sep = ""), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
      
      ## add gene symbol names
      exp.temp <- cbind(Symbol = anno.temp, exp.temp)
      exp.temp <- exp.temp <- exp.temp[exp.temp$Gene.symbol != "",]
      exp.temp <- aggregate(exp.temp[,-1], by = list(Symbol = exp.temp$Gene.symbol), FUN = mean, na.rm = TRUE)
      
      ## merge expression data
      if (i == 1) {
            expr <- exp.temp
      } else {
            expr <- merge(expr, exp.temp, by = "Symbol", all = TRUE)
      }
      
      ## print something
      cat("GSE18229-GPL", series[i], "data read, contains", ncol(exp.temp)-1, "samples\n")
      if (i == length(series)) {
            cat("All data read, in total contains", ncol(expr)-1, "samples\n")
      }
}
row.names(expr) <- expr$Symbol
expr$Symbol <- NULL
expr <- expr[,sort(colnames(expr))]
rm(anno.temp, exp.temp, i)

### read sample titles
for (i in 1:length(series)) {
      ## read sample titles
      sample.temp <- read.table(paste("data/GSE/GSE18229-GPL", series[i], ".txt", sep = ""), stringsAsFactors = FALSE, nrow = 2)
      sample.temp <- data.frame(t(sample.temp[,-1]))
      colnames(sample.temp) <- c("titles", "samples")
      
      ## merge sample titles
      if (i == 1) {
            sample <- sample.temp
      } else {
            sample <- rbind(sample, sample.temp)
      }
}
rm(sample.temp, i, series)

### read clinical data
clin <- read.table("data/GSE/pheno.GSE18229.txt", sep = "\t", header = TRUE, fill = TRUE)
clin <- clin[,c("GEO.array.names", "Overall.Survival.Event..0.alive..1.DOD.or.DOC.", "Overall.suvival.months", "Grade", "PAM50.calls.UNC337")]
colnames(clin) <- c("titles", "overall_survival_event", "overall_suvival_months", "Grade", "PAM50")
# filter out NA
idx <- apply(clin, 1, FUN = function(x) any(is.na(x)))
clin <- clin[!idx,]
## add sample name
clin <- merge(sample, clin, by = "titles")
row.names(clin) <- clin$samples; clin$titles <- NULL; clin$samples <- NULL
clin <- clin[sort(row.names(clin)),]

### filter out samples without clinical information
expr <- expr[,row.names(clin)]
rm(sample)

### save data
save.image("data/GSE/UNC_GSE18229_clean.RData")
