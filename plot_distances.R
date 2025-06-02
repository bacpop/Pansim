library(ggplot2)
library(reticulate)

path_to_python="/Users/shorsfield/miniforge3/envs/biopython/bin/python"
use_python(path_to_python)

py_run_string("
import numpy as np
from scipy.optimize import curve_fit

# fit asymptotic curve using exponential decay
def negative_exponential(x, b0, b1, b2): # based on https://isem-cueb-ztian.github.io/Intro-Econometrics-2017/handouts/lecture_notes/lecture_10/lecture_10.pdf and https://www.statforbiology.com/articles/usefulequations/
    return b0 * (1 - np.exp(-b1 * x)) + b2
    
# fit asymptotic curve using exponential decay
def negative_exponential_2(x, b0, b1, b2): # based on https://isem-cueb-ztian.github.io/Intro-Econometrics-2017/handouts/lecture_notes/lecture_10/lecture_10.pdf and https://www.statforbiology.com/articles/usefulequations/
    return b0 - (b0 - b1) * np.exp(-b2 * x)

def fit_curve(x, y):
  popt, pcov = curve_fit(negative_exponential, x, y, p0=[1.0, 1.0, 0.0], bounds=([0.0, 0.0, 0.0], [1.0, np.inf, 1.0]))
  b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))
  return popt.tolist(), [b0_err, b1_err, b2_err]
  
def fit_curve2(x, y):
  popt, pcov = curve_fit(negative_exponential_2, x, y, p0=[1.0, 0.0, 1.0], bounds=([0.0, 0.0, 0.0], [1.0, 1.0, np.inf]))
  b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))
  return popt.tolist(), [b0_err, b1_err, b2_err]
")

#plot pairwise distances
#filename <- "pan_mu_0.05_core_mu_0.019_prop_0.1_recomb_0.0_competition_true"
#in_dir <- "gridsearch/commit_5487550 [latest_core]/files/"
#out_dir <- "gridsearch/commit_5487550 [latest_core]/figures/"
in_dir <- "gridsearch/test2/"
out_dir <- "gridsearch/test2/"
filenames <- list.files(path = in_dir, pattern = "\\.tsv$")

for (filename in filenames)
{
  #filename="prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.05_rate_genes1_0.5_prop_genes2_0.1_threads_4_mu_core_0.019_1.tsv"
  df <- read.csv(paste(in_dir, filename, sep = ""), sep = "\t", header = FALSE)
  colnames(df) <- c("Core", "Accessory")
  
  adj_core = (-3/4) * log(1 - (4/3 * max(df$Core)))
  
  params <- py$fit_curve2(df$Core, df$Accessory)
  b0 <- params[[1]][[1]]
  b1 <- params[[1]][[2]]
  b2 <- params[[1]][[3]]
  b0_err <- params[[2]][[1]]
  b1_err <- params[[2]][[2]]
  b2_err <- params[[2]][[3]]

  p <- ggplot(df, aes(x = Core, y = Accessory)) + geom_point(alpha=0.2, colour="#0c589c") +
    labs(x = "Core", y = "Accessory") +
    geom_density_2d(aes(color = ..level..), bins = 15) +
    theme_light() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=18), legend.position ="none") +
    stat_function(fun = function(x) b0 - (b0 - b1) * exp(-b2 * x),
                  #fun = function(x) b0 * (1 - exp(-b1 * x)) + b2,
                  color = "red", linewidth = 1) +
      annotate("text", x = min(df$Core), y = max(df$Accessory), hjust = 0, vjust = 1,
               label = sprintf("Fitted: b0 = %.2f, b1 = %.2f, b2 = %.2f\nb0_err = %.2f, b1_err = %.2f, b2_err = %.2f", b0, b1, b2, b0_err, b1_err, b2_err),
               size = 5)
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

