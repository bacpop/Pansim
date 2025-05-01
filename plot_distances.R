library(ggplot2)

#plot pairwise distances
#filename <- "pan_mu_0.05_core_mu_0.019_prop_0.1_recomb_0.0_competition_true"
in_dir <- "gridsearch/commit_996b2e4/files/"
out_dir <- "gridsearch/commit_996b2e4/figures/"
filenames <- list.files(path = in_dir, pattern = "\\.tsv$")


for (filename in filenames)
{
  filename <- "pansim/distances.tsv"
  in_dir <- ""
  out_dir <- ""
  df <- read.csv(paste(in_dir, filename, sep = ""), sep = "\t", header = FALSE)
  colnames(df) <- c("Core", "Accessory")
  
  p <- ggplot(df, aes(x = Core, y = Accessory)) + geom_point(alpha=0.2, colour="#0c589c") +
    labs(x = "Core", y = "Accessory") +
    theme_light() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=18))
  p
  
  prefix <- sub(pattern = "(.*)\\..*$", replacement = "\\1", filename)
  ggsave(paste(out_dir, "ngen100_npop1000_", prefix, ".png", sep = ""))
}

# plot change in average pairwise distance
df <- read.csv("pansim/default_same_start_ngen_500_popsize_2000_per_gen.tsv", sep = "\t", header = FALSE)
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

ggsave("core_distance_per_gen_same_start_ngen500_popsize_2000_1_approx_avgfreq.png", plot=p.core, width=8, height=6)
ggsave("accessory_distance_per_gen_same_start_ngen500_popsize_2000_1_approx_avgfreq.png", plot=p.acc, width=8, height=6)
ggsave("core_distance_stddev_per_gen_same_start_ngen500_popsize_2000_1_approx_avgfreq.png", plot=p.core.std, width=8, height=6)
ggsave("accessory_distance_stddev_per_gen_same_ngen500_popsize_2000_1_approx_avgfreq.png", plot=p.acc.std, width=8, height=6)

