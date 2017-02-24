# =============================================== #
#              trimodal_survival.R                #
# =============================================== #
# corelate trimodal genes with survival data

source("code/select_source.R")
load("data/metabric/metabric_trimodal.maxExp.RData")
load("data/metabric/expr.max.RData")


# --------------------------- clean survival data ---------------------------- #
exp_discovery <- as.matrix(exp_discovery.max)
exp_validation <- as.matrix(exp_validation.max)
exp_normal <- as.matrix(exp_normal.max)
exp_discovery <- exp_discovery[,sort(colnames(exp_discovery))]
exp_validation <- exp_validation[,sort(colnames(exp_validation))]
TI_discovery_normal <- cbind(TI_discovery_normal.maxExp, data.frame(z_12 = NA, p_12 = NA, z_23 = NA, p_23 = NA, grade_12 = NA, grade_23 = NA, odd_12 = NA, odd_23 = NA))
TI_validation_normal <- cbind(TI_validation_normal.maxExp, data.frame(z_12 = NA, p_12 = NA, z_23 = NA, p_23 = NA, grade_12 = NA, grade_23 = NA, odd_12 = NA, odd_23 = NA))

### read survival data ###
anno_discovery <- read.table("data/metabric/table_S2_revised.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno_validation <- read.table("data/metabric/table_S3_revised.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

### clean data
anno_discovery$T <- anno_discovery$T/365.25
anno_validation$T <- anno_validation$T/365.25
anno_discovery$last_follow_up_status[anno_discovery$T > 20] <- 0
anno_validation$last_follow_up_status[anno_validation$T > 20] <- 0
anno_discovery$T[anno_discovery$T > 20] <- 20
anno_validation$T[anno_validation$T > 20] <- 20
# remove follow up status NA
# rm.idx_discovery <- (is.na(anno_discovery$last_follow_up_status)) | (anno_discovery$last_follow_up_status == "d-o.c.")
# rm.idx_validation <- (is.na(anno_validation$last_follow_up_status)) | (anno_validation$last_follow_up_status == "d-o.c.")
# anno_discovery <- anno_discovery[!rm.idx_discovery,]
# anno_validation <- anno_validation[!rm.idx_validation,]
# exp_discovery <- exp_discovery[,!rm.idx_discovery]
# exp_validation <- exp_validation[,!rm.idx_validation]
# factor follow up status
anno_discovery$last_follow_up_status[anno_discovery$last_follow_up_status %in% c("a", "d-o.c.")] <- 0
anno_validation$last_follow_up_status[anno_validation$last_follow_up_status %in% c("a", "d-o.c.")] <- 0
anno_discovery$last_follow_up_status[anno_discovery$last_follow_up_status %in% c("d", "d-d.s.")] <- 1
anno_validation$last_follow_up_status[anno_validation$last_follow_up_status %in% c("d", "d-d.s.")] <- 1
anno_discovery$last_follow_up_status <- as.numeric(anno_discovery$last_follow_up_status)
anno_validation$last_follow_up_status <- as.numeric(anno_validation$last_follow_up_status)


# -------------------------- correlate with survival data ------------------------------ #
pdf(file = "report/trimodal_survival.pdf", width = 8, height = 8)
par(mfrow = c(2,2))

ids_trimodal <- which(TI_discovery_normal$pi1 > 0.05 & TI_discovery_normal$pi3 > 0.05 & TI_validation_normal$pi1 > 0.05 & TI_validation_normal$pi3 > 0.05)
keep_id <- c()

# validate with survival and annotation information
for (i in ids_trimodal) {
      cat("Analyzing", i, "\n")
      
      cat_exp_discovery <- categorize(exp_discovery[i,], TI_discovery_normal$cutoff12[i], TI_discovery_normal$cutoff23[i])
      TI_discovery_normal[i,c("z_12", "p_12", "z_23", "p_23")] <- plot_surv(anno_discovery, cat_exp_discovery, title = rownames(exp_discovery)[i], other_covariate = NULL, plot = FALSE)
      TI_discovery_normal[i,c("grade_12", "grade_23", "odd_12", "odd_23")] <- check_anno(cat_exp_discovery, anno_discovery$grade, anno_discovery$Pam50Subtype == "Basal")
  
      cat_exp_validation <- categorize(exp_validation[i,],TI_validation_normal$cutoff12[i],TI_validation_normal$cutoff23[i])
      TI_validation_normal[i,c("z_12", "p_12", "z_23", "p_23")] <- plot_surv(anno_validation, cat_exp_validation, title = rownames(exp_validation)[i], other_covariate = NULL, plot = FALSE)
      TI_validation_normal[i,c("grade_12", "grade_23", "odd_12", "odd_23")] <- check_anno(cat_exp_validation, anno_validation$grade, anno_validation$Pam50Subtype == "Basal")
  
      pval_surv <- c(TI_discovery_normal[i,c("p_12","p_23")], TI_validation_normal[i,c("p_12","p_23")])
      pval_grade <- c(TI_discovery_normal[i,c("grade_12","grade_23")], TI_validation_normal[i,c("grade_12","grade_23")])
      
      ### look if survival and grade enrichment is
      if (all(pval_surv < 0.1) & all(pval_grade < 0.1)) {
            keep_id <- c(keep_id,i)
            cat(i, rownames(exp_discovery)[i], " is significant!\n", sep = "")
    
            trimodal <- trimodal_detect(exp_discovery[i,], exp_normal[i,], plot = TRUE, naive = FALSE, title = rownames(exp_discovery)[i])
            plot_surv(anno_discovery, cat_exp_discovery, title = "Discovery set", other_covariate = NULL, plot = TRUE)
    
            trimodal <- trimodal_detect(exp_validation[i,], exp_normal[i,], plot = TRUE, naive = FALSE, title = rownames(exp_validation)[i])
            plot_surv(anno_validation, cat_exp_validation, title = "Validation set", other_covariate = NULL, plot = TRUE)
      }
}
dev.off()

# -------------------------------- multivariate analysis ----------------------------------- #
i <- which("ORMDL3" == rownames(exp_normal))
covariates <- c("ERBB2")

anno_discovery$ERBB2 <- exp_discovery["ERBB2",]
cat_exp_discovery <- categorize(exp_discovery[i,], TI_discovery_normal$cutoff12[i], TI_discovery_normal$cutoff23[i])
tmp <- plot_surv(anno_discovery, cat_exp_discovery, rownames(exp_discovery)[i], covariates, plot = TRUE, print = TRUE)
coef <- coef(summary(lm(grade ~ ., data=data.frame(cat_exp=relevel(factor(cat_exp_discovery), ref="middle"), anno_discovery[,c("grade",covariates)]))))
coef[,4] <- coef[,4]/2
coef[coef[,1]<0,4] <- 1-coef[coef[,1]<0,4]
coef

anno_validation$ERBB2 <- exp_validation["ERBB2",]
anno_validation$Pam50Subtype[anno_validation$Pam50Subtype == "NC"] <- NA
cat_exp_validation <- categorize(exp_validation[i,],TI_validation_normal$cutoff12[i],TI_validation_normal$cutoff23[i])
tmp <- plot_surv(anno_validation, cat_exp_validation, rownames(exp_validation)[i], covariates, plot = TRUE, print = TRUE)
coef <- coef(summary(lm(grade ~ ., data=data.frame(cat_exp=relevel(factor(cat_exp_validation), ref="middle"), anno_validation[,c("grade",covariates)]))))
coef[,4] <- coef[,4]/2
coef[coef[,1]<0, 4] <- 1-coef[coef[,1]<0,4]
coef

anno_discovery$ESR1 <- exp_discovery["ESR1",]
anno_discovery$PGR <- exp_discovery["PGR",]
anno_discovery$ERBB2 <- exp_discovery["ERBB2",]
anno_discovery$Pam50Subtype[grepl("Lum",anno_discovery$Pam50Subtype)] <- "LumAB"
covariates <- c("ERBB2", "Pam50Subtype", "age_at_diagnosis", "stage", "ESR1", "PGR", "lymph_nodes_positive")

tmp <- plot_surv(anno_discovery, cat_exp_discovery, rownames(exp_discovery)[i], covariates, plot = TRUE, print = TRUE)

# -------------- save results ------------- #
TI_discovery_normal <- TI_discovery_normal[keep_id,]
TI_validation_normal <- TI_validation_normal[keep_id,]
rownames(TI_discovery_normal)
write.table(TI_discovery_normal, file="data/metabric/trimodal_survival_discovery.txt", quote = FALSE,row.names = TRUE, col.names = TRUE, sep="\t")
write.table(TI_validation_normal, file="data/metabric/trimodal_survival_validation.txt", quote = FALSE,row.names = TRUE, col.names = TRUE, sep="\t")
