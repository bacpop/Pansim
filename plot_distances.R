library(ggplot2)

#plot pairwise distances
df <- read.csv("distances.tsv", sep = "\t", header = FALSE)
colnames(df) <- c("Core", "Accessory")

p <- ggplot(df, aes(x = Core, y = Accessory)) + geom_point(alpha=0.2, colour="#0c589c") +
  labs(x = "Core", y = "Accessory") +
  theme_minimal()
p

# plot change in average pairwise distance
df <- read.csv("pansim/default_ngen_100_per_gen.tsv", sep = "\t", header = FALSE)
colnames(df) <- c("Avg_Core", "Std_Core", "Avg_Accessory", "Std_Accessory")
df$ID <- seq.int(nrow(df))
df$min_Core <- df$Avg_Core - df$Std_Core
df$max_Core <- df$Avg_Core + df$Std_Core
df$min_Acc <- df$Avg_Accessory - df$Std_Accessory
df$max_Acc <- df$Avg_Accessory + df$Std_Accessory

p.acc <- ggplot(df, aes(x = ID, y = Avg_Accessory)) + geom_line(colour="#0c589c") +
  geom_ribbon(aes(ymin = min_Acc, ymax = max_Acc), fill = "#0c589c", alpha = 0.2) +
  labs(x = "Generation", y = "Avg. Accessory Distance") +
  theme_light()
p.acc

p.acc.std <- ggplot(df) + geom_line(aes(x=ID, y=Std_Accessory), colour="#0c589c", linetype="dashed") +
  labs(x = "Generation", y = "StdDev. Accessory Distance") +
  theme_light()
p.acc.std

p.core <- ggplot(df, aes(x = ID, y = Avg_Core)) + geom_line(colour="#0c589c") +
  geom_ribbon(aes(ymin = min_Core, ymax = max_Core), fill = "#0c589c", alpha = 0.2) +
  labs(x = "Generation", y = "Avg. Core Distance") +
  theme_light()
p.core

p.core.std <- ggplot(df) + geom_line(aes(x=ID, y=Std_Core), colour="#0c589c", linetype="dashed") +
  labs(x = "Generation", y = "StdDev. Core Distance") +
  theme_light()
p.core.std

ggsave("core_distance_per_gen.png", plot=p.core, width=8, height=6)
ggsave("accessory_distance_per_gen.png", plot=p.acc, width=8, height=6)
ggsave("core_distance_stddev_per_gen.png", plot=p.core.std, width=8, height=6)
ggsave("accessory_distance_stddev_per_gen.png", plot=p.acc.std, width=8, height=6)
