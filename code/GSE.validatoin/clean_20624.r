# ============================================ #
#               clean_20624.R                  #
# ============================================ #

# -------------------- read data ------------------------ #
series <- c("885", "887", "1390", "5325")

### read expression data
for (i in 1:length(series)) {
      ## read expression and annotation data      
      exp.temp <- read.table(paste("data/GSE/GSE20624-GPL", series[i], "_series_matrix.txt", sep = ""), stringsAsFactors = FALSE, header = TRUE, comment.char = "!", row.names = 1, fill = TRUE)
      anno.temp <- read.table(paste("data/GSE/GPL", series[i], ".anno.txt", sep = ""), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
      
      ## add gene symbol names
      exp.temp <- cbind(Symbol = anno.temp, exp.temp)
      exp.temp <- exp.temp[exp.temp$Gene.symbol != "",]
      exp.temp <- aggregate(exp.temp[,-1], by = list(Symbol = exp.temp$Gene.symbol), FUN = mean, na.rm = TRUE)
      
      ## merge expression data
      if (i == 1) {
            expr <- exp.temp
      } else {
            expr <- merge(expr, exp.temp, by = "Symbol", all = TRUE)
      }
      
      ## print something
      cat("GSE20624-GPL", series[i], "data read, contains", ncol(exp.temp)-1, "samples\n")
      if (i == length(series)) {
            cat("All data read, in total contains", ncol(expr)-1, "samples\n")
      }
}
row.names(expr) <- expr$Symbol
expr$Symbol <- NULL
expr <- expr[,sort(colnames(expr))]
rm(anno.temp, exp.temp, i, series)

### read clinical data
clin <- read.csv("data/GSE/pheno.GSE20624.csv")
clin <- clin[clin$GEO.individual.microarray.Accession.number != "",]
row.names(clin) <- clin$GEO.individual.microarray.Accession.number
clin$GEO.individual.microarray.Accession.number <- NULL
clin <- clin[sort(row.names(clin)),]
# filter out NA
colnames(clin)[colnames(clin) == "PAM50.Call"] <- "PAM50"
clin <- clin[(!is.na(clin$Grade)) & (!is.na(clin$PAM50)),]

expr <- expr[,row.names(clin)]

### save data
save.image("data/GSE/UNC_GSE20624_clean.RData")
