# =========================================== #
#                 GSE20624.R                  #
# =========================================== #

source("code/select_source.R")

### ------------------------ read data -------------------------- ###
# potential hits
load("data/metabric/metabric_trimodal.maxExp.RData")
TI_metabric <- (TI_discovery_normal.maxExp + TI_validation_normal.maxExp) / 2
hits_normal <- rownames(TI_metabric)

# read breast cancer data
load("data/GSE/UNC_GSE20624_clean.RData")
colnames(clin)[colnames(clin) == "DOD_and_OS_months"] = "T"
colnames(clin)[colnames(clin) == "Dead_of_Disease"] = "last_follow_up_status"
clin$T = clin$T / 12

results <- matrix(-1, ncol = 9, nrow = length(hits_normal), dimnames = list(hits_normal, c("TI", "z_12", "p_12", "z_23", "p_23", "grade_12", "grade_23", "odd_12", "odd_23")))

#########  investigate survival data  ###############a

pdf("~/project/ORMDL3/metabric/trimodal/GSE20624.pdf", width = 6, height = 8)
par(mfrow = c(3, 2))

for (hit in hits_normal) {
      # get expression data
      if (!hit %in% rownames(expr)) {
            next
      }
      exp_hit = unlist(expr[hit, ])

      # plot trimodal distributino
      trimodal = list(TI = -1, cutoff12 = quantile(exp_hit, TI_metabric[hit, "pi1"], na.rm = T), cutoff23 = quantile(exp_hit, 1 - TI_metabric[hit, "pi3"], na.rm = T))

      # plot survival data
      results[hit, "TI"] = trimodal$TI
      cat_exp = categorize(exp_hit, trimodal$cutoff12, trimodal$cutoff23)
      results[hit, c("z_12", "p_12", "z_23", "p_23")] = plot_surv(clin, cat_exp, "GSE20624", other = NULL)

      # grade and ER grouping
      results[hit, c("grade_12", "grade_23", "odd_12", "odd_23")] <- check_anno(cat_exp, clin$Grade, clin$PAM50 == "Basal")
}

dev.off()

##########  multivariate analysis  ###############
hit <- "ORMDL3"
exp_hit <- unlist(expr[hit, ])
clin$ERBB2 <- unlist(expr["ERBB2", ])
covariates <- c("ERBB2")

# plot trimodal distribution
trimodal <- list(TI = -1, cutoff12 = quantile(exp_hit, TI_metabric[hit, "pi1"], na.rm = T), cutoff23 = quantile(exp_hit, 1 - TI_metabric[hit, "pi3"], na.rm = T))

# plot survival data
cat_exp <- categorize(exp_hit, trimodal$cutoff12, trimodal$cutoff23)
tmp <- plot_surv(clin, cat_exp, hit, other = covariates, print = T, plot = T)

# grade
coef = coef(summary(lm(Grade ~ ., data = data.frame(cat_exp = relevel(factor(cat_exp), ref = "middle"), clin[, c("Grade", covariates)]))))
coef[, 4] = coef[, 4] / 2
coef[coef[, 1] < 0, 4] = 1 - coef[coef[, 1] < 0, 4]
coef

###########  save results  ######################
# results = as.data.frame(results)
# write.table(results, file = "~/project/ORMDL3/metabric/trimodal/GSE20624.txt", quote = F, row.names = T, col.names = T, sep = "\t")
