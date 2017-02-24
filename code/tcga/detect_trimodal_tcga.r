# ================================================ #
#               detect_trimodal.R                  #
# ================================================ #
source("code/trimodal_detect.R")

# --------------------- read data -------------------------- #
load("data/tcga/unc.edu_BRCA_IlluminaHiSeq_RNASeqV2.geneExp.RData")

TI_tumor_normal <- as.data.frame(matrix(-1, ncol=11, nrow=dim(expr.tumor)[1], dimnames = list(rownames(expr.tumor), c("id", "TI", "mu1", "mu2", "mu3", "sigma", "pi1", "pi2", "pi3", "cutoff12", "cutoff23"))))

##########  detect trimodal distribution  ##############
error.lines <- numeric(0)
for (i in 1:dim(expr.tumor)[1]) {
      
      trimodal <- trimodal_detect(expr.tumor[i,], expr.normal[i,], naive = FALSE)
      
      if (!is.null(trimodal)) {
            if (is.list(trimodal)) {
                  TI_tumor_normal[i,] <- c(i,trimodal$TI, trimodal$mu, trimodal$sigma, trimodal$pi, trimodal$cutoff12, trimodal$cutoff23)
            } else if (trimodal == 0) {
                  error.lines <- c(error.lines, i)
                  cat(i, "contains errors\n")
            }
      }
      if (i %% 100 == 0) cat(i, "is done!\n")
}

############# save data ############################
save(TI_tumor_normal, file = "data/tcga/tcga_trimodal.RData")

