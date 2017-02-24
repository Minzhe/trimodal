# ====================================================================== #
#                           CNV.nearby.R                                 #
# ====================================================================== #
# Analyze the copy number variation of ERBB2 and ORMDL3.


### ---------------------- 1. discovery set ----------------------- ###

############ read cnv and cna data #############
cnv_discovery <- read.table("data/metabric/EGAD00010000214/discovery_CNV_CBS.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cna_discovery <- read.table("data/metabric/EGAD00010000213/discovery_CNA_CBS.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


###### cnv analysis #####

### look for ORMDL3 
cnv_discovery.ORMDL3 <- cnv_discovery[cnv_discovery$chrom == 17,]
cnv_discovery.ORMDL3 <- cnv_discovery.ORMDL3[cnv_discovery.ORMDL3$loc.end > 38000000 & cnv_discovery.ORMDL3$loc.start < 38200000,]
cnv_discovery.ORMDL3 <- cnv_discovery.ORMDL3[,c("METABRIC_ID", "seg.mean")]
cnv_discovery.ORMDL3 <- aggregate(cnv_discovery.ORMDL3$seg.mean, by = list(METABRIC_ID = cnv_discovery.ORMDL3$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cnv_discovery.ORMDL3)[2] <- "seg.mean.ORMDL3"


### look for ERBB2
cnv_discovery.ERBB2 <- cnv_discovery[cnv_discovery$chrom == 17,]
cnv_discovery.ERBB2 <- cnv_discovery.ERBB2[cnv_discovery.ERBB2$loc.end > 37800000 & cnv_discovery.ERBB2$loc.start < 38000000,]
cnv_discovery.ERBB2 <- cnv_discovery.ERBB2[,c("METABRIC_ID", "seg.mean")]
cnv_discovery.ERBB2 <- aggregate(cnv_discovery.ERBB2$seg.mean, by = list(METABRIC_ID = cnv_discovery.ERBB2$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cnv_discovery.ERBB2)[2] <- "seg.mean.ERBB2"

### merge cnv data
cnv_discovery.nearby <- merge(cnv_discovery.ORMDL3, cnv_discovery.ERBB2, by = "METABRIC_ID")

###### cna analysis ########

### look for ORMDL3 
cna_discovery.ORMDL3 <- cna_discovery[cna_discovery$chrom == 17,]
cna_discovery.ORMDL3 <- cna_discovery.ORMDL3[cna_discovery.ORMDL3$loc.end > 38000000 & cna_discovery.ORMDL3$loc.start < 38200000,]
cna_discovery.ORMDL3 <- cna_discovery.ORMDL3[,c("METABRIC_ID", "seg.mean")]
cna_discovery.ORMDL3 <- aggregate(cna_discovery.ORMDL3$seg.mean, by = list(METABRIC_ID = cna_discovery.ORMDL3$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cna_discovery.ORMDL3)[2] <- "seg.mean.ORMDL3"


### look for ERBB2
cna_discovery.ERBB2 <- cna_discovery[cna_discovery$chrom == 17,]
cna_discovery.ERBB2 <- cna_discovery.ERBB2[cna_discovery.ERBB2$loc.end > 37800000 & cna_discovery.ERBB2$loc.start < 38000000,]
cna_discovery.ERBB2 <- cna_discovery.ERBB2[,c("METABRIC_ID", "seg.mean")]
cna_discovery.ERBB2 <- aggregate(cna_discovery.ERBB2$seg.mean, by = list(METABRIC_ID = cna_discovery.ERBB2$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cna_discovery.ERBB2)[2] <- "seg.mean.ERBB2"

### merge cna data
cna_discovery.nearby <- merge(cna_discovery.ORMDL3, cna_discovery.ERBB2, by = "METABRIC_ID")
rm(cna_discovery.ORMDL3, cna_discovery.ERBB2)
rm(cnv_discovery.ORMDL3, cnv_discovery.ERBB2)


### -------------------------- 2. validation set ------------------------------- ###

########### read data ##############
cnv_validation <- read.table("data/metabric/EGAD00010000216/validation_CNV_CBS.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cna_validation <- read.table("data/metabric/EGAD00010000215/validation_CNA_CBS.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


########### cnv analysis ############

### look for ORMDL3 
cnv_validation.ORMDL3 <- cnv_validation[cnv_validation$chrom == 17,]
cnv_validation.ORMDL3 <- cnv_validation.ORMDL3[cnv_validation.ORMDL3$loc.end > 38000000 & cnv_validation.ORMDL3$loc.start < 38200000,]
cnv_validation.ORMDL3 <- cnv_validation.ORMDL3[,c("METABRIC_ID", "seg.mean")]
cnv_validation.ORMDL3 <- aggregate(cnv_validation.ORMDL3$seg.mean, by = list(METABRIC_ID = cnv_validation.ORMDL3$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cnv_validation.ORMDL3)[2] <- "seg.mean.ORMDL3"


### look for ERBB2
cnv_validation.ERBB2 <- cnv_validation[cnv_validation$chrom == 17,]
cnv_validation.ERBB2 <- cnv_validation.ERBB2[cnv_validation.ERBB2$loc.end > 37800000 & cnv_validation.ERBB2$loc.start < 38000000,]
cnv_validation.ERBB2 <- cnv_validation.ERBB2[,c("METABRIC_ID", "seg.mean")]
cnv_validation.ERBB2 <- aggregate(cnv_validation.ERBB2$seg.mean, by = list(METABRIC_ID = cnv_validation.ERBB2$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cnv_validation.ERBB2)[2] <- "seg.mean.ERBB2"

### merge cnv data
cnv_validation.nearby <- merge(cnv_validation.ORMDL3, cnv_validation.ERBB2, by = "METABRIC_ID")


############ cna analysis ###############

### look for ORMDL3 
cna_validation.ORMDL3 <- cna_validation[cna_validation$chrom == 17,]
cna_validation.ORMDL3 <- cna_validation.ORMDL3[cna_validation.ORMDL3$loc.end > 38000000 & cna_validation.ORMDL3$loc.start < 38200000,]
cna_validation.ORMDL3 <- cna_validation.ORMDL3[,c("METABRIC_ID", "seg.mean")]
cna_validation.ORMDL3 <- aggregate(cna_validation.ORMDL3$seg.mean, by = list(METABRIC_ID = cna_validation.ORMDL3$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cna_validation.ORMDL3)[2] <- "seg.mean.ORMDL3"


### look for ERBB2
cna_validation.ERBB2 <- cna_validation[cna_validation$chrom == 17,]
cna_validation.ERBB2 <- cna_validation.ERBB2[cna_validation.ERBB2$loc.end > 37800000 & cna_validation.ERBB2$loc.start < 38000000,]
cna_validation.ERBB2 <- cna_validation.ERBB2[,c("METABRIC_ID", "seg.mean")]
cna_validation.ERBB2 <- aggregate(cna_validation.ERBB2$seg.mean, by = list(METABRIC_ID = cna_validation.ERBB2$METABRIC_ID), FUN = function(x) x[which.max(abs(x))])
colnames(cna_validation.ERBB2)[2] <- "seg.mean.ERBB2"

### merge cna data
cna_validation.nearby <- merge(cna_validation.ORMDL3, cna_validation.ERBB2, by = "METABRIC_ID")
rm(cnv_validation.ORMDL3, cnv_validation.ERBB2)
rm(cna_validation.ORMDL3, cna_validation.ERBB2)

############## plot something ####################
pdf("~/project/ORMDL3/metabric/trimodal/figures/cnvcna.pdf", width=6, height=8)
par(mfrow = c(3,2))

plot(x = cnv_discovery.nearby$seg.mean.ORMDL3, y = cnv_discovery.nearby$seg.mean.ERBB2, main = "Discovery set CNV data", xlab = "ORMDL3", ylab = "ERBB2", pch = 19, col = "red", cex = 0.8, xlim = c(-0.6, 1.5), ylim = c(-0.6, 1.5))
points(x = cna_discovery.nearby$seg.mean.ORMDL3, y = cna_discovery.nearby$seg.mean.ERBB2, pch = 19, col = "red", cex = 0.8)
plot(x = cnv_validation.nearby$seg.mean.ORMDL3, y = cnv_validation.nearby$seg.mean.ERBB2, main = "validation set CNV data", xlab = "ORMDL3", ylab = "ERBB2", pch = 19, col = "red", cex = 0.8, xlim = c(-0.6, 1.2), ylim = c(-0.6, 1.2))
points(x = cna_validation.nearby$seg.mean.ORMDL3, y = cna_validation.nearby$seg.mean.ERBB2, pch = 19, col = "red", cex = 0.8)

dev.off()
