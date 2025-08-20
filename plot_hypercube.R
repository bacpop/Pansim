library(ggplot2)
library(patchwork)

#plot hypercube parameter correlations
filename <- ""
df <- read.csv(filename, sep = "\t", header = TRUE)

#params <- c("b0", "b1", "b2", "b0_err", "b1_err", "b2_err", "mean_acc")
params <- c("b0", "b1", "b1_err", "mean_acc")
pairs <- combn(params, 2, simplify = FALSE)

plot_list <- list()

i <- 1
for (i in seq_along(pairs))
{
  x_col <- pairs[[i]][1]
  y_col <- pairs[[i]][2]
  
  corr <- cor(df[[x_col]], df[[y_col]], use = "complete.obs")
  title <- paste(x_col, "vs", y_col, "- Cor:", round(corr, 2))
  
  p <- ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(alpha = 0.6) +
    labs(title = title, x = x_col, y = y_col) +
    theme_light()
  p
  ggsave(paste("param_correlation_", x_col, "_" , y_col,".png", sep=""), plot=p, width=8, height=6)
  #plot_list[[i]] <- p
}

# Combine and display all plots in a grid
p.all <- wrap_plots(plotlist = plot_list)

ggsave("hypercube_parameter_correlation.png", plot=p.all, width=8, height=6)

