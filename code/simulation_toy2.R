# Install the following packages:
#   install_github("stephenslab/flashr@v0.5-6", build_vignettes = TRUE)
#   install.packages("PMA")
#   install.packages("ssvd")
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("impute", version = "3.8")


library(flashr) # Empirical bayes matrix factorization
library(PMA) # Penalized matrix decomposition + Sparse PCA
library(impute)
library(reshape2)
library(ggplot2)
source("code/helpers.R")
###################################################################################
##   5.1.1. Single factor example with sparse loadings and a non-sparse factor   ##
###################################################################################
# specify the dimensions
n <- 200
p <- 300
k <- 1

# initialize other parameters
pi_0 <- 0.9 # loading sparsity options: 0.9, 0.3, 0
tau <- 1 # corresponsing noise precision options: 1, 1/16, 1/25

#sparse loadings from a mixture Gaussian model
sigma_squared <- c(0.25, 0.5, 1, 2, 4)
num_gm <- length(sigma_squared)
mus <- rep(0, length(sigma_squared))
mixture_prob <- rep(1/5, num_gm)

num_sim <- 50
flash_error <- rep(NA, num_sim)
flash_v2_error <- rep(NA, num_sim)
flash_backfit_error <- rep(NA, num_sim)
PMD_error <- rep(NA, num_sim)
PMD2_error <- rep(NA, num_sim)
SSVD_error <- rep(NA, num_sim)
SVD_error <- rep(NA, num_sim)

for (i in 1:num_sim) {
  print("start")
  print(i)
  
  set.seed(i)
  mixture_ids <- sample(1:length(sigma_squared), n, replace = TRUE, prob = mixture_prob)
  sparsity_ids <- sample(c(0, 1), prob = c(pi_0, 1 - pi_0), size = n, replace = TRUE)
  LL <- sparsity_ids * rnorm(n = n, mean = mus[mixture_ids], sd = sqrt(sigma_squared[mixture_ids]))
  
  # non-sparse factor
  FF <- rnorm(p) # non-sparse
  
  # original matrix: Y = LF^T + E 
  LF <- LL %*% t(FF)
  E <- matrix(rnorm(n * p, mean = 0, sd = sqrt(1/tau)), nrow = n)
  Y <- LF + E
  
  # run flash version 1
  data <- flash_set_data(Y)
  flash_out <- flash(data, Kmax = 1, verbose=TRUE, ebnm_fn = ebnm_ash)
  
  # run flash version 2
  flash_v2_out <- flash(data, Kmax = 1, verbose=TRUE, ebnm_fn = ebnm_pn)
  
  # run backfitting
  flash_backfit_out = flash(Y, f_init = flash_v2_out, backfit = TRUE, greedy = FALSE, verbose = TRUE)
  
  # run PMD + CV
  PMD_cv_out <- PMD.cv(Y, type = "standard", sumabss = seq(0.1, 0.9, len=20))
  PMD_out <- PMD(Y, type = "standard", sumabs = PMD_cv_out$bestsumabs, K = 1, v = PMD_cv_out$v.init)
  
  
  PMD_cv_out2 <- PMD.cv(Y, type = "ordered", sumabsus = seq(1, sqrt(n), len = 10))
  PMD_out2 <- PMD(Y, type ="ordered", sumabsu=PMD_cv_out2$bestsumabsu, K = 1, v = PMD_cv_out2$v.init)
  
  # run Sparse PCA
  # SPC_cv_out <- SPC.cv(Y)
  # out.orth <- SPC(Y, sumabsv = SPC_cv_out$bestsumabsv, K = 1, orth=TRUE)
  
  # run SSVD
  SSVD_out <- ssvd(Y, method = "method", r = 1, non.orth=TRUE)
  
  # run SVD
  SVD_out <- svd(Y, nu = 1, nv = 1)
  
  # compute RMSE
  flash_error[i] <- sqrt(mean((LF- flash_get_lf(flash_out))^2))
  flash_v2_error[i] <- sqrt(mean((LF- flash_get_lf(flash_v2_out))^2))
  flash_backfit_error[i] <- sqrt(mean((LF- flash_get_lf(flash_backfit_out))^2))
  PMD_error[i] <- sqrt(mean((LF - PMD_out$d * PMD_out$u  %*% t(PMD_out$v))^2))
  PMD2_error[i] <- sqrt(mean((LF - PMD_out2$d * PMD_out2$u  %*% t(PMD_out2$v))^2))
  SSVD_error[i] <- sqrt(mean((LF - SSVD_out$d * SSVD_out$u  %*% t(SSVD_out$v))^2))
  SVD_error[i] <- sqrt(mean((LF - SVD_out$d * SVD_out$u  %*% t(SVD_out$v))^2))
}

df <- data.frame(#flash = flash_error,
  EBMF_PN = flash_v2_error,
  EBMF_backfit = flash_backfit_error,
  PMD = PMD_error,
  #                 PMD2 = PMD2_error,
  SSVD = SSVD_error)#,
#SVD = SVD_error)
long <- melt(df)
names(long) <- c("Method", "value")

#write.csv(long, "output/csv_files/toy2.csv")

toy2 <- ggplot(data = long) + 
  geom_boxplot(aes(y = value, color = Method)) + 
  scale_y_sqrt() + ylab("RRMSE") + 
  mytheme

ggsave(filename = "toy2.png", path = "output/figures", plot = toy2, height = 5, width = 5)



# potentially useful plots
figure_LF <- plot_matrix(LF)
ggsave(filename = "figure_LF2.png", path = "output/figures", plot = figure_LF, height = 5, width = 5)

figure_LF_flash <- plot_matrix(flash_get_lf(flash_v2_out))
ggsave(filename = "figure_LF_flash2.png", path = "output/figures", plot = figure_LF_flash, height = 5, width = 5)

figure_flash_reconstruct_error <- plot_matrix(LF - flash_get_lf(flash_v2_out))
ggsave(filename = "figure_flash_reconstruct_error2.png", path = "output/figures", plot = figure_flash_reconstruct_error, height = 5, width = 5)

figure_PMD_reconstruct_error <- plot_matrix(LF - PMD_out$d * PMD_out$u  %*% t(PMD_out$v))
ggsave(filename = "figure_PMD_reconstruct_error2.png", path = "output/figures", plot = figure_PMD_reconstruct_error, height = 5, width = 5)

figure_SSVD_reconstruct_error <- plot_matrix(LF - SSVD_out$d * SSVD_out$u  %*% t(SSVD_out$v))
ggsave(filename = "figure_SSVD_reconstruct_error2.png", path = "output/figures", plot = figure_SSVD_reconstruct_error, height = 5, width = 5)







