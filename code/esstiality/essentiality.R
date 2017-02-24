############### load data  #############

hit <- "ORMDL3"
nearby <- "ERBB2"

load("data/CCLE/CCLE_Expression_clean.RData")
load("data/GARP/garp.score.RData")
rm(expr.ccle.raw, garp.data.raw)

############### clean data ##############
### ccle expression data
expr.ccle <- expr.ccle[,grep("BREAST", colnames(expr.ccle))]
colnames(expr.ccle) <- gsub("_BREAST", "", colnames(expr.ccle))
expr.ccle <- expr.ccle[,sort(colnames(expr.ccle))]

### grap score data
colnames(garp.score) <- toupper(gsub("[.]", "", colnames(garp.score)))

### garp cancer cancer cell lines
garp.cellLines$Cell.line <- toupper(gsub("-| ", "", garp.cellLines$Cell.line))
row.names(garp.cellLines) <- garp.cellLines$Cell.line; garp.cellLines$Cell.line <- NULL

### select breast cancer cell lines
brca.cell <- row.names(garp.cellLines)[garp.cellLines$Tumor.type == "Breast"]
brca.cell <- intersect(brca.cell, colnames(expr.ccle))

expr.ccle <- as.matrix(expr.ccle[,brca.cell])
garp.score <- as.matrix(garp.score[,brca.cell])
garp.cellLines <- garp.cellLines[brca.cell,]

############ gene essentiality and expression  #############

### all subtypes
pdf("~/project/ORMDL3/metabric/trimodal/figures/gene.essentiality.pdf", width=5, height=5)
par(mfrow = c(2,2))

plot_garp <- function(gene, correct = T) {
      x <- expr.ccle[gene,]
      y <- garp.score[gene,]
      if (correct == T) {
            x <- x - predict(lm(x~factor(garp.cellLines$Subtype)))
            y <- y - predict(lm(y~factor(garp.cellLines$Subtype)))
      }
      plot(x, y, main = gene, xlab = "expression", ylab = "Gene essentiality", pch=19, col = "gold")
      lines(x[order(x)], predict(lm(y~x))[order(x)], lwd = 3, col = "purple")
      text(min(x)*0.95 + max(x)*0.05, min(y)*0.95 + max(y)*0.05, paste("Pear. corr.=", round(cor(x,y), d=2)), col = "purple", pos = 4)
}

plot_garp(nearby)
plot_garp(hit)

plot_garp("MYC")
plot_garp("CCND1")
plot_garp("FAM83B")

plot_garp("TP53")
plot_garp("BRCA1")
plot_garp("BRCA2")
plot_garp("PTEN")
plot_garp("ATM")
plot_garp("RAD50")

dev.off()

### Luminal

pdf("~/project/ORMDL3/metabric/trimodal/figures/gene.essentiality.Luminal.pdf", width=5, height=5)
par(mfrow = c(2,2))

plot_garp_subtype <- function(gene, subtype) {
      x <- expr.ccle[gene, row.names(garp.cellLines)[garp.cellLines$Subtype == subtype]]
      y <- garp.score[gene, row.names(garp.cellLines)[garp.cellLines$Subtype == subtype]]

      plot(x, y, main = gene, xlab = paste("expression", subtype), ylab = "Gene essentiality", pch=19, col = "gold")
      lines(x[order(x)], predict(lm(y~x))[order(x)], lwd = 3, col = "purple")
      text(min(x)*0.95 + max(x)*0.05, min(y)*0.95 + max(y)*0.05, paste("Pear. corr.=", round(cor(x,y), d=2)), col = "purple", pos = 4)
}

plot_garp_subtype(nearby, "Luminal")
plot_garp_subtype(hit, "Luminal")

plot_garp_subtype("MYC", "Luminal")
plot_garp_subtype("CCND1", "Luminal")
plot_garp_subtype("FAM83B", "Luminal")

plot_garp_subtype("TP53", "Luminal")
plot_garp_subtype("BRCA1", "Luminal")
plot_garp_subtype("BRCA2", "Luminal")
plot_garp_subtype("PTEN", "Luminal")
plot_garp_subtype("ATM", "Luminal")
plot_garp_subtype("RAD50", "Luminal")

dev.off()


### BasalA/B
garp.cellLines$Subtype[garp.cellLines$Subtype %in% c("Basal A", "Basal B")] <- "BasalA/B"
pdf("~/project/ORMDL3/metabric/trimodal/figures/gene.essentiality.Basal.pdf", width=5, height=5)
par(mfrow = c(2,2))

plot_garp_subtype(nearby, "BasalA/B")
plot_garp_subtype(hit, "BasalA/B")

plot_garp_subtype("MYC", "BasalA/B")
plot_garp_subtype("CCND1", "BasalA/B")
plot_garp_subtype("FAM83B", "BasalA/B")

plot_garp_subtype("TP53", "BasalA/B")
plot_garp_subtype("BRCA1", "BasalA/B")
plot_garp_subtype("BRCA2", "BasalA/B")
plot_garp_subtype("PTEN", "BasalA/B")
plot_garp_subtype("ATM", "BasalA/B")
plot_garp_subtype("RAD50", "BasalA/B")

dev.off()

### HER2
garp.cellLines$Subtype[garp.cellLines$Subtype %in% c("Basal A", "Basal B")] <- "HER2"
pdf("~/project/ORMDL3/metabric/trimodal/figures/gene.essentiality.HER2.pdf", width=5, height=5)
par(mfrow = c(2,2))

plot_garp_subtype(nearby, "HER2")
plot_garp_subtype(hit, "HER2")

plot_garp_subtype("MYC", "HER2")
plot_garp_subtype("CCND1", "HER2")
plot_garp_subtype("FAM83B", "HER2")

plot_garp_subtype("TP53", "HER2")
plot_garp_subtype("BRCA1", "HER2")
plot_garp_subtype("BRCA2", "HER2")
plot_garp_subtype("PTEN", "HER2")
plot_garp_subtype("ATM", "HER2")
plot_garp_subtype("RAD50", "HER2")

dev.off()

