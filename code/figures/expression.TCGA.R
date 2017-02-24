# ====================================================================== #
#                          expression.TCGA.R                             #
# ====================================================================== #
# Analyze the gene expression of ERBB2 and ORMDL3.

################ read TCGA data ####################
expr.raw <- read.table("data/tcga/unc.edu_BRCA_IlluminaHiSeq_RNASeqV2.geneExp.tsv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")

################ clean date ######################
expr <- round(log(expr.raw + 1), 6)
type <- sapply(strsplit(colnames(expr), ".", fixed = TRUE), "[[", 4)
tumor_id <- grep("^0", type)
normal_id <- grep("^1", type)
expr.tumor <- expr[,tumor_id]
expr.normal <- expr[,normal_id]

################ ORMDL3 and ERBB2 ##################
expr.tumor.nearby <- expr.tumor[grep("ERBB2|ORMDL3", row.names(expr.tumor)),][-1,]
expr.normal.nearby <- expr.normal[grep("ERBB2|ORMDL3", row.names(expr.normal)),][-1,]

############# make plots #############
pdf("~/project/ORMDL3/metabric/trimodal/figures/TCGA.expression.pdf", width = 6, height = 8)
par(mfrow = c(3,2))

plot(as.numeric(expr.tumor.nearby["ORMDL3",]), as.numeric(expr.tumor.nearby["ERBB2",]), pch = 19, col = "red", cex = 0.8, main = "TCGA expression data", xlab = "ORMDL3", ylab = "ERBB2")
points(as.numeric(expr.normal.nearby["ORMDL3",]), as.numeric(expr.normal.nearby["ERBB2",]), pch = 19, col = "blue", cex = 0.8)

dev.off()
