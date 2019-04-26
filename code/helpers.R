# helper functions
#\bibitem{plot_matrix} https://rpubs.com/lgadar/matrix-visualizations


plot_matrix <- function(X) {
  long <- melt(X)
  long<-long[long$value!=0, ]
  
  p <- ggplot(long, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient2()
}