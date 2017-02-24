# =========================================== #
#              select_source.R                #
# =========================================== #
library(survival)
source("code/trimodal_detect.R")

### ------------------------ survival curve --------------------------- ###
# correlate trimodal genes with survival data
plot_surv <- function(anno, cat_exp, title, other_covariate = NULL, plot = TRUE, print = FALSE) { 
  
      cat_exp_factor <- factor(cat_exp, labels = 1:3,levels = c("middle", "low", "high"))
      fit <- survfit(Surv(anno$T, anno$last_follow_up_status) ~ cat_exp_factor) 
      
      form <- "Surv(anno$T, anno$last_follow_up_status) ~ relevel(cat_exp_factor, ref = 1)"
      if (!is.null(other_covariate)) {
            form <- paste(c(form, paste("+anno$", other_covariate, sep = "")), collapse = "")
      }
      
      cox <- summary(coxph(as.formula(form)))$coef
      cox[,5] <- cox[,5] / 2
      cox[cox[,4] < 0, 5] <- 1 - cox[cox[,4] < 0, 5]
      colnames(cox)[5] <- "pval"
      
      if (print == TRUE) {
            cox <- cbind(cox,matrix(0, ncol = 2, nrow = dim(cox)[1]))
            colnames(cox)[c(2,6:7)] <- c("HR","L.HR","U.HR")
            cox[,6] <- cox[,2] / exp(1.96*cox[,3])
            cox[,7] <- cox[,2] * exp(1.96*cox[,3])
            print(cox)
      }
      cox12 <- cox["relevel(cat_exp_factor, ref = 1)2", c("z","pval")]
      cox23 <- cox["relevel(cat_exp_factor, ref = 1)3", c("z","pval")]

  
      if (plot == TRUE) {
            tmp <- as.numeric(unique(cat_exp_factor))
            plot(fit, col = tmp[order(tmp)], mark = 19, main = paste("Kaplan-Meier plot for", title), xlab = "Time (years)", ylab = "Survival probability")
            text(0.01 * max(anno$T, na.rm = TRUE), 0.05, paste("p12=", pretty_num(cox12[2])), col = "red", pos = 4)
            text(0.01 * max(anno$T, na.rm = TRUE), 0.15, paste("p23=", pretty_num(cox23[2])), col = "green", pos = 4)
            text(0.01 * max(anno$T, na.rm = TRUE), 0.25, "Lower group", col = "red", pos = 4)
            text(0.01 * max(anno$T, na.rm = TRUE), 0.35, "Middle group", col = "black", pos = 4)
            text(0.01 * max(anno$T, na.rm = TRUE), 0.45, "Higher group", col = "green", pos = 4)
      }
  
      # p values are two-way
      c(cox12, cox23)
}


# ------------------------------ check enrichment ---------------------------------- #
# check enrichment of high and low fraction in high grade tumor and Triple-negative/basal subtype
# grade should be a numeric vector of 1, 2, 3 as well as other elements
# tn should be a boolean vector with T equal to triple-negative or basal subtype
check_anno <- function(cat_exp, grade, tn) {
      
      grade[is.na(grade) | (!grade %in% c(1,2,3))] <- 4
      count <- as.matrix(table(cat_exp,factor(grade,labels=1:4,levels=1:4)))

      if (!"low" %in% cat_exp) {
            pt_low_middle <- NA
            or_low_middle <- NA
      } else {
            pt_low_middle <- prop.trend.test(count["low",1:3], count["low",1:3] + count["middle",1:3], c(1,2,3))$p.value/2
            
            if (cor(count["low",1:3]/(count["low",1:3]+count["middle",1:3]),1:3,use="complete") < 0) {
                  pt_low_middle <- 1 - pt_low_middle/2
            }
            or_low_middle <- oddsratio(sum(cat_exp == "low" & tn,na.rm = TRUE), sum(cat_exp=="middle" & tn,na.rm = TRUE), sum(cat_exp == "low" & (!tn), na.rm = TRUE), sum(cat_exp=="middle" & (!tn), na.rm = TRUE))
      }
      
      if (!"high" %in% cat_exp) {
            pt_high_middle <- NA
            or_high_middle <- NA
      } else {
            pt_high_middle <- prop.trend.test(count["high",1:3], count["high",1:3] + count["middle",1:3], c(1,2,3))$p.value/2
            
            if (cor(count["high",1:3]/(count["high",1:3]+count["middle",1:3]),1:3,use="complete")<0) {
                  pt_high_middle <- 1 - pt_high_middle/2
            }
            or_high_middle <- oddsratio(sum(cat_exp == "high" & tn,na.rm = TRUE), sum(cat_exp == "middle" & tn,na.rm = TRUE), sum(cat_exp == "high" & (!tn),na.rm = TRUE), sum(cat_exp == "middle" & (!tn), na.rm = TRUE))
      }
      
      c(pt_low_middle, pt_high_middle, or_low_middle, or_high_middle)
}

### ----------------------- auxiliary function ------------------------- ###
categorize <- function(x, cutoff1, cutoff2) {
      tmp <- x
      tmp[!is.na(tmp)] <- "middle"
      tmp[x<cutoff1] <- "low"
      tmp[x>cutoff2] <- "high"
      tmp
}

pretty_num <- function(x) {
      pwr <- ceiling(-log10(x))
      paste(round(10^pwr*x, d=3), "e-", pwr, sep = "")
}

oddsratio<-function(x11,x12,x21,x22) {
      c <- 0.5
      x11 <- x11 + c
      x12 <- x12 + c
      x21 <- x21 + c
      x22 <- x22 + c
      or <- log(x22*x11/x21/x12)
      se <- sqrt(1/x11+1/x12+1/x21+1/x22)
      pnorm(-or/se)
}