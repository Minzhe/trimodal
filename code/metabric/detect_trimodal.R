# ================================================ #
#               detect_trimodal.R                  #
# ================================================ #
# trimodal detection of metabric data set.

setwd("/qbrc/home/mzhang/projects/trimodal")
source("code/trimodal_detect.R")

# -------------- detect trimodal for max expression data ---------------- #
##### discovery set #####
### read data ###
load("data/metabric/expr.max.RData")
exp_discovery <- exp_discovery.max
exp_validation <- exp_validation.max
exp_normal <- exp_normal.max

TI_discovery_normal <- as.data.frame(matrix(-1, ncol = 11, nrow = dim(exp_discovery)[1], dimnames = list(rownames(exp_discovery), c("id", "TI", "mu1", "mu2", "mu3", "sigma", "pi1", "pi2", "pi3", "cutoff12", "cutoff23"))))

### detect trimodal distribution
for (i in 1:dim(exp_discovery)[1]) {
      trimodal <- trimodal_detect(exp_discovery[i,], exp_normal[i,], naive = FALSE)
      
      if (!is.null(trimodal)) {
            if(is.list(trimodal)) {
                  TI_discovery_normal[i,] <- c(i, trimodal$TI, trimodal$mu, trimodal$sigma, trimodal$pi, trimodal$cutoff12, trimodal$cutoff23)
            } else if (trimodal == 0) {
                  error.lines <- c(error.lines, i)
                  cat(i, "contains errors\n")
            }
      }
      if (i %% 100 == 0) cat(i, "is done!\n")
}

TI_discovery_normal.maxExp <- TI_discovery_normal

# -------------- detect trimodal for validation set ---------------- #
##### validation set #####
TI_validation_normal <- as.data.frame(matrix(-1, ncol = 11, nrow = dim(exp_validation)[1], dimnames = list(rownames(exp_validation), c("id", "TI", "mu1", "mu2", "mu3", "sigma", "pi1", "pi2", "pi3", "cutoff12", "cutoff23"))))

### detect trimodal distribution ###
for (i in 1:dim(exp_validation)[1]) {
      trimodal <- trimodal_detect(exp_validation[i,], exp_normal[i,], naive = FALSE)
      
      if (!is.null(trimodal)) {
            if(is.list(trimodal)) {
                  TI_validation_normal[i,] <- c(i, trimodal$TI, trimodal$mu, trimodal$sigma, trimodal$pi, trimodal$cutoff12, trimodal$cutoff23)
            } else if (trimodal == 0) {
                  error.lines <- c(error.lines, i)
                  cat(i, "contains errors\n")
            }
      }
      if (i %% 100 == 0) cat(i, "is done!\n")
}

TI_validation_normal.maxExp <- TI_validation_normal


# -------------------------- save data ----------------------------- #
save(TI_discovery_normal.maxExp, TI_validation_normal.maxExp, file = "data/metabric/metabric_trimodal.maxExp.RData")
