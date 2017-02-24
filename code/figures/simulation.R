# =========================================== #
#               simulation.R                  #
# =========================================== #
# Simulation of trimodal distribution

setwd("/qbrc/home/mzhang/projects/trimodal")
source("code/select_source.R")

pdf(file = "report/simulation.pdf", 5, 5)
x <- c(rnorm(100,-4,1), rnorm(250,0,1), rnorm(150,3,1))
x0 <- rnorm(50,0,1)
trimodal_detect(x, x0, plot = TRUE, title = "simulated example")
dev.off()