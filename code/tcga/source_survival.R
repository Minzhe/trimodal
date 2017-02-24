# get expression data
load("~/iproject/virus/data/expression/tcga/expression.RData")
subtype = "BRCA"
exp_subtype = exp[[subtype]]

# read survival data
surv_files = list.files(
      "~/iproject/trimodal/data/tcga/survival/",
      pattern = "childrens",
      full.names = T
)
surv_data = read.csv(surv_files[grepl(subtype, surv_files)], sep = "\t", stringsAsFactors =
                           F)
tmp = table(surv_data$tissue_source_site)
tmp = names(tmp[tmp > 20])
surv_data = surv_data[surv_data$tissue_source_site %in% tmp, ]
surv_data$tissue_source_site = as.factor(surv_data$tissue_source_site)

surv_data$last_contact_days_to[is.na(surv_data$last_contact_days_to)] =
      surv_data$death_days_to[is.na(surv_data$last_contact_days_to)]
surv_data = surv_data[surv_data$vital_status != "" &
                            (!is.na(surv_data$last_contact_days_to)), ]
surv_data$vital_status = surv_data$vital_status == "Dead"
surv_data$bcr_patient_barcode = gsub("-", "\\.", surv_data$bcr_patient_barcode)
rownames(surv_data) = paste(surv_data$bcr_patient_barcode, ".01", sep =
                                  "")
surv_data = surv_data[rownames(surv_data) %in% rownames(exp_subtype), ]

colnames(surv_data)[colnames(surv_data) == "last_contact_days_to"] = "T"
colnames(surv_data)[colnames(surv_data) == "vital_status"] = "last_follow_up_status"

surv_data$TN = F
surv_data$TN[surv_data$er_status_by_ihc == "Negative" &
                   surv_data$her2_status_by_ihc == "Negative" &
                   surv_data$pr_status_by_ihc == "Negative"] = T
indeterminate = c("", "Equivocal", "Indeterminate")
surv_data$TN[surv_data$er_status_by_ihc %in% indeterminate |
                   surv_data$her2_status_by_ihc %in% indeterminate |
                   surv_data$pr_status_by_ihc %in% indeterminate] = NA

load("~/iproject/trimodal/data/tcga/survival/cleaned_BECN1_TCGA_1007.RData")
surv_data$grade = NA
surv_data$grade[surv_data$sample %in% newd$PID[newd$pathReport_Grade == "I" &
                                                     (!is.na(newd$pathReport_Grade))]] = 1
surv_data$grade[surv_data$sample %in% newd$PID[newd$pathReport_Grade == "II" &
                                                     (!is.na(newd$pathReport_Grade))]] = 2
surv_data$grade[surv_data$sample %in% newd$PID[newd$pathReport_Grade == "III" &
                                                     (!is.na(newd$pathReport_Grade))]] = 3

#load("~/iproject/virus/data/mutation/tcga/mutation.RData")
#mut_subtype=mut[[subtype]]
