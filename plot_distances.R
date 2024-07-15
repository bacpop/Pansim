library(ggplot2)

df <- read.csv("distances.tsv", sep = "\t", header = FALSE)
colnames(df) <- c("Core", "Accessory")

p <- ggplot(df, aes(x = Core, y = Accessory)) + geom_point(alpha=0.2, colour="#0c589c") +
  labs(x = "Core", y = "Accessory") +
  theme_minimal()
p
