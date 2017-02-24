# ============================================ #
#               clean_20685.R                  #
# ============================================ #

library(Biobase)

# -------------------- read data ------------------------ #
load("data/GSE/BRCA_raw327TW_GSE20685.RData")


# ------------------- clean data -------------------- #
expr <- GSE20685.processed

clin <- pData(GSE20685[[1]])
names(clin)[names(clin) == "characteristics_ch1.3"] <- "last_follow_up_status"
names(clin)[names(clin) == "characteristics_ch1.4"] <- "T"
clin <- clin[,c("last_follow_up_status", "T")]

clin$last_follow_up_status <- as.character(clin$last_follow_up_status)
clin$last_follow_up_status <- gsub("event_death: ", "", clin$last_follow_up_status)
clin$last_follow_up_status <- as.numeric(clin$last_follow_up_status)

clin$T <- as.character(clin$T)
clin$T <- gsub("follow_up_duration \\(years\\): ", "", clin$T)
clin$T <- as.numeric(clin$T)


# ------------------- save data ----------------------- #
rm(GSE20685.processed, GSE.genes, GSE20685, cb.data, kmc, res, cb.TCGA)
save.image("data/GSE/UNC_GSE20685_clean.RData")
