############# read data #############
garp.data.raw <- read.table("~/project/ORMDL3/GARP/GARP-score.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
garp.cellLines <- read.table("~/project/ORMDL3/GARP/GARP.celllines.txt", skip = 2, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote = "")

############ clean data #############
garp.score <- garp.data.raw[,-c(1,3,4)]
garp.score <- aggregate(garp.score[,-1], by = list(Gene.name = garp.score$Gene.name), FUN = mean)
row.names(garp.score) <- garp.score$Gene.name
garp.score$Gene.name <- NULL

garp.cellLines <- garp.cellLines[1:72,]


########### save data ###############
save.image("~/project/ORMDL3/GARP/garp.score.RData")
