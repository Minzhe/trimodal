################ read data ######################
expr.ccle.raw <- read.csv("project/ORMDL3/CCLE/CCLE_Expression_Entrez_2012-10-18.csv")


############### clean data ######################
expr.ccle <- expr.ccle.raw[-(1:2),-grep("^X.", colnames(expr.ccle))]
expr.ccle <- expr.ccle[,-2]
expr.ccle$Description <- gsub("expression ", "", expr.ccle$Description)
expr.ccle <- expr.ccle[expr.ccle$Description != "",]
expr.ccle$Accession <- NULL
expr.ccle <- expr.ccle[!duplicated(expr.ccle$Description),]
row.names(expr.ccle) <- expr.ccle$Description
expr.ccle$Description <- NULL


############## save data #######################
save.image("~/project/ORMDL3/CCLE/CCLE_Expression_clean.RData")

