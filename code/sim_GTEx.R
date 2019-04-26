library(ashr)
library(flashr)
dat <- readRDS("data/MatrixEQTLSumStats.Portable.Z.rds")
zscore <- dat$strong.z
Y <- t(zscore)
data <- flash_set_data(Y)
f_greedy <- flash_add_greedy(data,Kmax=60)
f_greedy_bf <- flash_backfit(data,f_greedy)
#f_greedy = flash_add_greedy(data,Kmax=60,var_type = "by_column",ash_param=list(method = "fdr"))
#f_greedy_bf = flash_backfit(data,f_greedy,var_type = "by_column",ash_param=list(method = "fdr"))
saveRDS(f_greedy, file = "output/gflashvarcol.rds")
saveRDS(f_greedy_bf, file = "output/bflashvarcol.rds")

b_flash = readRDS("output/gflashvarcol.rds")
#load("output/GTExdata/gtexEQTL_zscore.rds")
ssY = sum(zscore^2)
K = dim(b_flash$EL)[2] -1
pve = (sapply(seq(1,K),function(x){ sum(b_flash$EL[,x]^2 %*% t(b_flash$EF[,x]^2)) }))/ssY
pve = pmax(round(pve,3),0.001)
dat = read.table('data/GTExColors.txt', sep = '\t', comment.char = '')
colordata = dat[c(1:6,9:18,21:23,26:30,32,33,35,36,38:53),1:2]
L = b_flash$EL[,1:14]
library(reshape2)
data_L = melt(L)
colnames(data_L) = c("tissue","loading","value")
library(ggplot2)
tissue_color = as.character(colordata[,2])
data_L$tissue = factor(data_L$tissue,levels = 1:44 ,labels = as.character(colordata[,1]) )
data_L$loading = factor(data_L$loading,levels = 1:14 ,labels = paste("Factor",1:14,"; pve:", pve[1:14]))
ggplot(data_L,aes(x = tissue,y = value,fill = factor(tissue) )) +
  geom_bar(stat = "identity",width = 0.6) +
  scale_fill_manual(values=tissue_color) +
  scale_x_discrete(labels = NULL) +
  theme_grey()+
  theme(legend.position="right", legend.text=element_text(size=9), axis.text.y = element_text(size = 5)) + 
  labs(title = "GTEx data", y = "factor values" ,x = "tissues", fill="tissue") +
  facet_wrap(~loading, ncol = 2, scales = "free_y") +
  guides(fill = guide_legend(ncol = 1, keyheight = 0.8, keywidth = 0.3))
ggsave("flashrGTEx1.png", path = "output/figures", width = 8, height = 11)
# the 27th factor is zero
L = b_flash$EL[,15:26]
library(reshape2)
data_L = melt(L)
colnames(data_L) = c("tissue","loading","value")
library(ggplot2)
tissue_color = as.character(colordata[,2])
data_L$tissue = factor(data_L$tissue,levels = 1:44 ,labels = as.character(colordata[,1]) )
data_L$loading = factor(data_L$loading,levels = 1:12 ,labels = paste("Factor",15:26,"; pve:", pve[15:26]))
ggplot(data_L,aes(x = tissue,y = value,fill = factor(tissue) )) +
  geom_bar(stat = "identity",width = 0.6) +
  scale_fill_manual(values=tissue_color) +
  scale_x_discrete(labels = NULL) +
  theme_grey()+
  theme(legend.position="right", legend.text=element_text(size=9), axis.text.y = element_text(size = 5)) + 
  labs(title = "GTEx data", y = "factor values" ,x = "tissues", fill="tissue") +
  facet_wrap(~loading, ncol = 2, scales = "free_y") +
  guides(fill = guide_legend(ncol = 1, keyheight = 0.8, keywidth = 0.3))
ggsave("flashrGTEx2.pdf", path = "output/figures", width = 8, height = 10)
