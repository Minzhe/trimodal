# ====================================================================== #
#                       expression.metabric.R                            #
# ====================================================================== #
# Analyze the gene expression of ERBB2 and ORMDL3.

################ read data ####################
load("data/metabric/expr.max.RData")

################ ORMDL3 and ERBB2 ##############
exp_discovery.nearby <- exp_discovery.max[c("ORMDL3", "ERBB2"),]
exp_validation.nearby <- exp_validation.max[c("ORMDL3", "ERBB2"),]
exp_normal.nearby <- exp_normal.max[c("ORMDL3", "ERBB2"),]

################ make plot ###################
pdf("~/project/ORMDL3/metabric/trimodal/expression.nearby.pdf", width = 6, height = 8)
par(mfrow = c(3,2))

plot(as.numeric(exp_discovery.nearby["ORMDL3",]), as.numeric(exp_discovery.nearby["ERBB2",]), pch = 19, col = "red", cex = 0.8, main = "Discovery set expression data", xlab = "ORMDL3", ylab = "ERBB2")
points(as.numeric(exp_normal.nearby["ORMDL3",]), as.numeric(exp_normal.nearby["ERBB2",]), pch = 19, col = "blue", cex = 0.8)

plot(as.numeric(exp_validation.nearby["ORMDL3",]), as.numeric(exp_validation.nearby["ERBB2",]), pch = 19, col = "red", cex = 0.8, main = "Validation set expression data", xlab = "ORMDL3", ylab = "ERBB2")
points(as.numeric(exp_normal.nearby["ORMDL3",]), as.numeric(exp_normal.nearby["ERBB2",]), pch = 19, col = "blue", cex = 0.8)

dev.off()