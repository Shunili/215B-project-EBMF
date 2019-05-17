# helper functions
#\bibitem{plot_matrix} https://rpubs.com/lgadar/matrix-visualizations

mytheme <- theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #plot.title = element_text(family = "serif", size = 16, hjust = 0.5), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        #axis.title.x = element_text(family = "serif", size = 15),
        axis.text.x = element_blank(),
        #axis.text.y = element_text(size = 14),
        axis.ticks.x=element_blank())

matrixtheme <- theme_classic() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        #plot.title = element_text(family = "serif", size = 16, hjust = 0.5), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())



plot_matrix <- function(X, title = " ") {
  long <- melt(X)
  long<-long[long$value!=0, ]
  
  p <- ggplot(long, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + ylab("samples (n)") + xlab("variables (p)")+ ggtitle(title)+
    scale_fill_gradient2() + matrixtheme
}
