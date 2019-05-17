library(flashr) # Empirical bayes matrix factorization
library(PMA) # Penalized matrix decomposition + Sparse PCA
library(ssvd)
library(impute)
library(reshape2)
library(ggplot2)
source("code/helpers.R")

###################################################################################
##   5.1.2. Sparse bi-cluster rank 3 example                                     ##
###################################################################################
# specify the dimensions

n <- 150
p <- 240
k <- 3
tau <- 1/4

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
  
  LL <- matrix(0, nrow = n, ncol = k)
  LL[1:10, 1] <- rnorm(10, mean = 0, sd = 2)
  LL[11:60, 2] <- rnorm(50, mean = 0, sd = 1)
  LL[61:150, 3] <- rnorm(90, mean = 0, sd = 0.5)

  
  FF <- matrix(0, nrow = p, ncol = k)
  FF[1:80, 1] <- rnorm(80, mean = 0, sd = 0.5)
  FF[81:160, 2] <- rnorm(80, mean = 0, sd = 1)
  FF[161:240, 3] <- rnorm(80, mean = 0, sd = 2)

  # original matrix: Y = LF^T + E 
  LF <- LL %*% t(FF)
  E <- matrix(rnorm(n * p, mean = 0, sd = sqrt(1/tau)), nrow = n)
  Y <- LF + E
  
  # run flash version 1
  data <- flash_set_data(Y)
  flash_out <- flash(data, Kmax = 3, verbose=TRUE, ebnm_fn = ebnm_ash)
  
  # run flash version 2
  flash_v2_out <- flash(data, Kmax = 1, verbose=TRUE, ebnm_fn = ebnm_pn)
  
  # run backfitting
  flash_backfit_out = flash(Y, f_init = flash_v2_out, backfit = TRUE, greedy = FALSE, verbose = TRUE)
  
  # run PMD + CV
  PMD_cv_out <- PMD.cv(Y, type = "standard", sumabss = seq(0.1, 0.9, len=20))
  PMD_out <- PMD(Y, type = "standard", sumabs = PMD_cv_out$bestsumabs, K = k, v = PMD_cv_out$v.init)
  
  
  PMD_cv_out2 <- PMD.cv(Y, type = "ordered", sumabsus = seq(1, sqrt(n), len = 10))
  PMD_out2 <- PMD(Y, type ="ordered", sumabsu=PMD_cv_out2$bestsumabsu, K = k, v = PMD_cv_out2$v.init)
  
  # run Sparse PCA
  # SPC_cv_out <- SPC.cv(Y)
  # out.orth <- SPC(Y, sumabsv = SPC_cv_out$bestsumabsv, K = 1, orth=TRUE)
  
  # run SSVD
  SSVD_out <- ssvd(Y, method = "method", r = 3)#, non.orth=TRUE)
  
  # run SVD
  SVD_out <- svd(Y, nu = 3, nv = 3)
  
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

bi_cluster <- ggplot(data = long) + 
  geom_boxplot(aes(y = value, color = Method)) + 
  scale_y_sqrt() + 
  ylab("RRMSE") + mytheme
ggsave(filename = "bi_cluster.png", path = "output/figures", plot =bi_cluster, height = 5, width = 5)

pp1<-plot_matrix(LF, title = "EBMF_PN")
pp2<-plot_matrix(flash_get_lf(flash_v2_out), title = "EBMF_PN")
pp3<-plot_matrix(abs(flash_get_lf(flash_backfit_out)), title = "EBMF_backfit")
pp4<-plot_matrix(abs(PMD_out$u  %*% t(PMD_out$v)), title = "PMD")
pp5<-plot_matrix(abs(SSVD_out$u  %*% t(SSVD_out$v)), title = "SSVD")

ggsave(filename = "pp1_bi_cluster.png", path = "output/figures", plot =pp1, height = 5, width = 5)
ggsave(filename = "pp2_bi_cluster.png", path = "output/figures", plot =pp1, height = 5, width = 5)
ggsave(filename = "pp3_bi_cluster.png", path = "output/figures", plot =pp3, height = 5, width = 5)
ggsave(filename = "pp4_bi_cluster.png", path = "output/figures", plot =pp4, height = 5, width = 5)
ggsave(filename = "pp5_bi_cluster.png", path = "output/figures", plot =pp5, height = 5, width = 5)

