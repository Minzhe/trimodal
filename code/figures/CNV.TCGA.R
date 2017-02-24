# =============================================== #
#                    CNV.TCGA.R                   #
# =============================================== #
# copy number variation analysis of ORMDL3 and ERBB2 in TCGA

################ read data ####################

cnv.data <- read.table("data/tcga/broad.mit.edu_BRCA_Genome_Wide_SNP_6.hg19.cna.seg", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


################ cnv analysis ##################
### ORMDL3
cnv.ORMDL3 <- cnv.data[cnv.data$Chromosome == 17,]
cnv.ORMDL3 <- cnv.ORMDL3[cnv.ORMDL3$End > 38077296 & cnv.ORMDL3$Start < 38083884,]
cnv.ORMDL3 <- cnv.ORMDL3[,c("Sample", "Segment_Mean")]
cnv.ORMDL3 <- aggregate(cnv.ORMDL3$Segment_Mean, by = list(Sample = cnv.ORMDL3$Sample), FUN = mean)
colnames(cnv.ORMDL3)[2] <- "Segment_Mean.ORMDL3"


### ERBB2
cnv.ERBB2 <- cnv.data[cnv.data$Chromosome == 17,]
cnv.ERBB2 <- cnv.ERBB2[cnv.ERBB2$End > 37856254 & cnv.ERBB2$Start < 37884915,]
cnv.ERBB2 <- cnv.ERBB2[,c("Sample", "Segment_Mean")]
cnv.ERBB2 <- aggregate(cnv.ERBB2$Segment_Mean, by = list(Sample = cnv.ERBB2$Sample), FUN = mean)
colnames(cnv.ERBB2)[2] <- "Segment_Mean.ERBB2"

### merge cnv data
cnv.nearby <- merge(cnv.ORMDL3, cnv.ERBB2, by = "Sample")

############ make plot ################
pdf("~/project/ORMDL3/metabric/trimodal/figures/TCGA.CNV.pdf", height = 8, width = 6)
par(mfrow = c(3,2))

plot(x = cnv.nearby$Segment_Mean.ORMDL3, y = cnv.nearby$Segment_Mean.ERBB2, main = "TCGA CNV data", xlab = "ORMDL3", ylab = "ERBB2", pch = 19, col = "red", cex = 0.8, xlim = c(-1, 3.8), ylim = c(-1, 3.8))

dev.off()
