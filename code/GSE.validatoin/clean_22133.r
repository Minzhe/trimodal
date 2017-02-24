# ============================================ #
#               clean_22133.R                  #
# ============================================ #

# -------------------- read data ------------------------ #
### read expression data
expr <- read.table("data/GSE/GSE22133-GPL5345_series_matrix.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, comment.char = "!", row.names = 1, fill = TRUE)
anno <- read.table("data/GSE/GSE22133.anno.txt", sep = "\t", header = TRUE, row.names = 1)
expr <- cbind(Symbol = anno$GSE22133_geneSymbol, expr)
expr <- aggregate(expr[,-1], by = list(Symbol = expr$Symbol), FUN = mean, na.rm = TRUE)
row.names(expr) <- expr$Symbol
expr$Symbol <- NULL
expr <- expr[,sort(colnames(expr))]
rm(anno)

### read clinical data
sample <- read.table("data/GSE/GSE22133-GPL5345_series_matrix.txt", sep = "\t", skip = 102, nrow = 1, stringsAsFactors = FALSE)
OS <- read.table("data/GSE/GSE22133-GPL5345_series_matrix.txt", sep = "\t", skip = 66, nrow = 1, stringsAsFactors = FALSE)
Death <- read.table("data/GSE/GSE22133-GPL5345_series_matrix.txt", sep = "\t", skip = 65, nrow = 1, stringsAsFactors = FALSE)
Grade <- read.table("data/GSE/GSE22133-GPL5345_series_matrix.txt", sep = "\t", skip = 69, nrow = 1, stringsAsFactors = FALSE)
Subtype <- read.table("data/GSE/GSE22133-GPL5345_series_matrix.txt", sep = "\t", skip = 72, nrow = 1, stringsAsFactors = FALSE)
clin <- data.frame(t(rbind(OS[-1], Death[-1], Grade[-1], Subtype[-1])), row.names = sample[,-1], stringsAsFactors = FALSE)
colnames(clin) <- c("OS", "Death", "Grade", "Subtype")
# clean data
clin$OS <- gsub("os: ", "", clin$OS)
clin$Death <- gsub("osbin: ", "", clin$Death)
clin$Grade <- gsub("grade: ", "", clin$Grade)
clin$Subtype <- gsub("pam50 classification:", "", clin$Subtype)
clin$OS <- as.numeric(clin$OS)
clin$Death <- as.numeric(clin$Death)
clin$Grade <- as.numeric(clin$Grade)
# filter out NA
clin <- clin[!is.na(clin$Grade),]
expr <- expr[, colnames(expr) %in% row.names(clin)]
rm(Death, OS, sample, Grade, Subtype)

### save data
save.image("data/GSE/UNC_GSE22133_clean.RData")
