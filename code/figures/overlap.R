# ================================================ #
#               detect_trimodal.R                  #
# ================================================ #

# -------------------- Metabric ----------------------- #

load("data/metabric/metabric_trimodal.maxExp.RData")
TI_discovery_normal <- TI_discovery_normal.maxExp
TI_validation_normal <- TI_validation_normal.maxExp

a <- sum(TI_discovery_normal[, 1] > -1 & TI_validation_normal[, 1] == -1)
b <- sum(TI_discovery_normal[, 1] == -1 & TI_validation_normal[, 1] > -1)
c <- sum(TI_discovery_normal[, 1] > -1 & TI_validation_normal[, 1] > -1)
d <- sum(TI_discovery_normal[, 1] == -1 & TI_validation_normal[, 1] == -1)

a
b
c
d
phyper(c, b + c, a + d, a + c, lower = FALSE, log.p = TRUE)

all_metabric <- rownames(TI_discovery_normal)
trimodal_metabric <- rownames(TI_discovery_normal)[TI_discovery_normal[, 1] > -1 & TI_validation_normal[, 1] > -1]

# --------------------- TCGA ----------------------- #
load("data/tcga/tcga_trimodal.RData")
geneNames <- sapply(strsplit(row.names(TI_tumor_normal), split = "|", fixed = TRUE), "[[", 1)
uni.idx <- match(unique(geneNames), geneNames)
TI_tumor_normal <- TI_tumor_normal[uni.idx,]
row.names(TI_tumor_normal) <- geneNames[uni.idx]
TI_tumor_normal <- TI_tumor_normal[row.names(TI_tumor_normal) != "?",]
all_tcga <- row.names(TI_tumor_normal)
trimodal_tcga <- row.names(TI_tumor_normal)[TI_tumor_normal[, 1] > -1]

# ------------- compare shared genes between metabric and tcga -------------- #
all_metabric <- all_metabric[all_metabric %in% all_tcga]
trimodal_metabric <- trimodal_metabric[trimodal_metabric %in% all_tcga]
all_tcga <- all_tcga[all_tcga %in% all_metabric]
trimodal_tcga <- trimodal_tcga[trimodal_tcga %in% all_metabric]

c <- sum(trimodal_metabric %in% trimodal_tcga)
a <- sum(!trimodal_metabric %in% trimodal_tcga)
b <- sum(!trimodal_tcga %in% trimodal_metabric)
d <- sum(!all_tcga %in% c(trimodal_metabric, trimodal_tcga))

a
b
c
d
phyper(c, b + c, a + d, a + c, lower = FALSE, log.p = TRUE)
