# ================================================== #
#                  clean_metabric.R                  #
# ================================================== #


# ------------------- read and clean data -------------------------- #
exp_discovery <- read.table("data/metabric/EGAD00010000210/discovery_ExpressionMatrix.txt", sep = "\t")
exp_validation <- read.table("data/metabric/EGAD00010000211/validation_ExpressionMatrix.txt", sep = "\t")
exp_normal <- read.table("data/metabric/EGAD00010000212/normals_ExpressionMatrix.txt", sep = "\t")
save(exp_discovery, exp_validation, exp_normal, file = "data/metabric/rawdata.RData")

gene.anno <- read.csv("data/metabric/illuminaHT-12.v3.anno.csv", row.names = 1, stringsAsFactors = FALSE)

### 1. filter out unigene ###
exp_discovery <- exp_discovery[row.names(exp_discovery) %in% row.names(gene.anno),]
exp_validation <- exp_validation[row.names(exp_validation) %in% row.names(gene.anno),]
exp_normal <- exp_normal[row.names(exp_normal) %in% row.names(gene.anno),]

### 2. reorder row.names ###
exp_discovery <- exp_discovery[sort(row.names(exp_discovery)),]
exp_validation <- exp_validation[sort(row.names(exp_validation)),]
exp_normal <- exp_normal[sort(row.names(exp_normal)),]

### 3. add Symbol ###
Symbol <- gene.anno[row.names(exp_discovery), "Symbol"]
exp_discovery <- cbind(Symbol = Symbol, exp_discovery)
exp_validation <- cbind(Symbol = Symbol, exp_validation)
exp_normal <- cbind(Symbol = Symbol, exp_normal)

### 4. collapse probes to gene symbol ###
### using mean
exp_discovery.mean <- aggregate(exp_discovery[,-1], by = list(factor(exp_discovery$Symbol)), FUN = mean)
exp_validation.mean <- aggregate(exp_validation[,-1], by = list(factor(exp_validation$Symbol)), FUN = mean)
exp_normal.mean <- aggregate(exp_normal[,-1], by = list(factor(exp_normal$Symbol)), FUN = mean)

row.names(exp_discovery.mean) <- exp_discovery.mean$Group.1
row.names(exp_validation.mean) <- exp_validation.mean$Group.1
row.names(exp_normal.mean) <- exp_normal.mean$Group.1
exp_discovery.mean$Group.1 <- NULL
exp_validation.mean$Group.1 <- NULL
exp_normal.mean$Group.1 <- NULL

### using max
exp_discovery.max <- aggregate(exp_discovery[,-1], by = list(factor(exp_discovery$Symbol)), FUN = max)
exp_validation.max <- aggregate(exp_validation[,-1], by = list(factor(exp_validation$Symbol)), FUN = max)
exp_normal.max <- aggregate(exp_normal[,-1], by = list(factor(exp_normal$Symbol)), FUN = max)

row.names(exp_discovery.max) <- exp_discovery.max$Group.1
row.names(exp_validation.max) <- exp_validation.max$Group.1
row.names(exp_normal.max) <- exp_normal.max$Group.1
exp_discovery.max$Group.1 <- NULL
exp_validation.max$Group.1 <- NULL
exp_normal.max$Group.1 <- NULL

save(exp_discovery.mean, exp_validation.mean, exp_normal.mean, file = "data/metabric/expr.mean.RData")
save(exp_discovery.max, exp_validation.max, exp_normal.max, file = "data/metabric/expr.max.RData")
