# =========================================== #
#            trimodal_detect.R                #
# =========================================== #

# Function to deteact trimodal distribution genes.

### ------------------- parametes-------------------- ###
# x: the expression vector of tumor samples
# x0: the expression vector of normal samples. x0 may not exist
# n.cvg: the number of LL values to look back for convergence
# pval.cvg: the pvalue cutoff for checking convergence
# max.iterate: the maximum number of iterates allowed
# plot: a bool variable to determine whether a density plot should be drawn
# naive: a bool variable to decide whether a naive algorithm should be used
# q_cutoff: the quantile cutoff for naive algorithm

### ----------------- return ------------------- ###
# TI:
# mu: vector of the mean expression of low, normal and high group
# sigma: vector of the standard deviation of the expression of low, normal and high group
# pi: vector of proportion of patients belong to low, normal and high group
# cutoff12: expression cutoff for low and normal group
# cutoff23: expression cutoff for normal and high group

### --------------------------- main function ---------------------------- ###
trimodal_detect <- function(x, x0, n.cvg = 10, pval.cvg = 0.1, max.iterate = 200, plot = FALSE, naive = FALSE, q_cutoff = 0.1, title="") {
  
      # clean data and get initial estimate
      x <- x[!is.na(x)]
      x0 <- x0[!is.na(x0)]
      
      sigma <- sd(x0)
      mu <- c(min(x), median(x0), max(x))
      pi <- c(0.1, 0.8, 0.1)
      LL <- c()
  
      ### EM algorithm
      if (naive == FALSE) { 
            ### start iteration
            repeat {
                  # calculate LL
                  LL <- c(LL, cal_LL(x, x0, mu, sigma, pi))
      
                  if (any(is.na(LL)) | any(abs(LL) == Inf)) return(NULL)
                  if (length(LL) > max.iterate) break
                  
                  # check whether LL has reached convergence
                  if (length(LL) > 2 * n.cvg) {

                        newer_LL <- tail(LL, n.cvg) 
                        older_LL <- head(tail(LL, 2*n.cvg), n.cvg)
                        
                        # no significant increase detected
                        if (sd(newer_LL) < 1e-10) break
                        if (t.test(newer_LL, older_LL, alternative = "greater")$p.val > pval.cvg) break
                  }
      
                  # get responsiblity
                  gamma <- sapply(1:3, function(i) pi[i]*dnorm_safe(x, mu[i], sigma))
                  gamma <- gamma / apply(gamma, 1, sum)
    
                  # update parameters
                  pi <- apply(gamma, 2, mean)
      
                  mu[1] <- sum(x*gamma[,1]) / sum(gamma[,1])
                  mu[2] <- (sum(x*gamma[,2]) + sum(x0)) / (sum(gamma[,2]) + length(x0))
                  mu[3] <- sum(x*gamma[,3]) / sum(gamma[,3])
                  if (any(apply(gamma,2,sum) == 0)) {return(0)}
                  
                  if (mu[1] > mu[2] & mu[2] <= mu[3]) {
                        mu[1] <- mu[2] <- (sum(gamma[,1:2]*x) + sum(x0)) / (sum(gamma[,1:2]) + length(x0))
                  } else if (mu[1] <= mu[2] & mu[2] > mu[3]) {
                        mu[2] <- mu[3] <- (sum(gamma[,2:3]*x) + sum(x0)) / (sum(gamma[,2:3]) + length(x0))
                  } else if (mu[1] > mu[2] & mu[2] > mu[3]) {
                        mu[1] <- mu[2] <- mu[3] <- (sum(gamma*x) + sum(x0)) / (sum(gamma) + length(x0))
                  }
  
                  tmp <- sapply(1:3, function(i) (x-mu[i])^2)
                  sigma <- sqrt((sum(tmp*gamma) + sum((x0-mu[2])^2)) / (sum(gamma) + length(x0)))
            }
    
            ### pi is the mixture proportion from EM 
            if (mu[1] == mu[2] || mu[2] == mu[3] || any(pi <= 0.01)) return(NULL)
            cutoff12 <- (mu[1]^2-mu[2]^2-log(pi[1]/pi[2])*2*sigma^2)/2/(mu[1]-mu[2])
            cutoff23 <- (mu[2]^2-mu[3]^2-log(pi[2]/pi[3])*2*sigma^2)/2/(mu[2]-mu[3])
            if (cutoff12 >= mu[2] | cutoff12 <= mu[1]) cutoff12 <- quantile(x, q_cutoff)
            if (cutoff23 <= mu[2] | cutoff23 >= mu[3]) cutoff23 <- quantile(x, 1-q_cutoff)
    
            ### pi is the empirical proportion calculated from cutoffs
            pi <-  c(sum(x < cutoff12), sum(x >= cutoff12 & x <= cutoff23), sum(x > cutoff23)) / length(x)
            if (any(pi <= 0.01)) return(NULL)
            TI <- min(mu[3] - mu[2], mu[2] - mu[1]) / sigma * min(pi)^0.5
            
      } else {
            cutoff12 <- quantile(x, q_cutoff)
            cutoff23 <- quantile(x, 1-q_cutoff)
            pi <- c(sum(x < cutoff12), sum(x >= cutoff12 & x <= cutoff23), sum(x > cutoff23)) / length(x) 
            TI <- min(pi[1], pi[3])
      }
  
      if (plot == TRUE) {
            ds_x <- density(x)
            ds_x0 <- density(x0)
            ds_x0$y <- ds_x0$y * length(x0) / length(x)
            x_span <- range(c(x,x0))
            x_span <- c(x_span[1] - (x_span[2] - x_span[1]) / 5, x_span[2] + (x_span[2] - x_span[1]) / 5) 
            yhei <- max(c(ds_x$y, ds_x0$y))
            
            plot(ds_x, ylim = c(0,yhei), xlim = x_span, xlab = "Expression level", lwd = 3, col = "red", main = paste("Density plot of expression levels\n", title, sep = ""))
            lines(ds_x0, lwd = 3, col = "blue")
            
            abline(v = cutoff12, lwd = 3,col = "purple", lty = 3)
            abline(v = cutoff23, lwd = 3,col = "purple", lty = 3)
            
            text(x_span[1], yhei*0.9, "Normal", col = "blue", pos = 4)
            text(x_span[1], yhei*0.8, "Tumor", col = "red", pos = 4)
            text(x_span[1], yhei*0.7, "Cutoff", col = "purple", pos = 4)
      }
  
      list(TI = TI, mu = mu, sigma = sigma, pi = pi, cutoff12 = cutoff12, cutoff23 = cutoff23)
}

### ----------------------- auxiliary function ------------------------- ###

dnorm_safe <- function(v, mean, sd) {
      v[(v-mean)/sd > 30] <- mean + sd * 30
      v[(v-mean)/sd < (-30)] <- mean - sd * 30
      dnorm(v,mean,sd)
}

cal_LL <- function(x, x0, mu, sigma, pi) {
      sum(log(sapply(1:3, function(j) dnorm_safe(x, mu[j], sigma)) %*% pi)) + sum(log(dnorm_safe(x0, mu[2], sigma)))
}