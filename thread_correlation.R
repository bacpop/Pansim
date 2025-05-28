library(dplyr)
library(tidyr)
library(ggplot2)

df_raw <- read.csv("threading_test.txt", sep = "\t", header = FALSE)

# Step 1: Split by group and assign iteration per group
df_split <- df_raw %>%
  group_by(V1) %>%
  mutate(iteration = row_number()) %>%
  ungroup()

# Step 2: Pivot to wide format — one column per group
df_wide <- df_split %>%
  pivot_wider(names_from = V1, values_from = V2)

# Step 3: Pivot back to long format for plotting
df_long <- df_wide %>%
  pivot_longer(cols = -iteration, names_to = "group", values_to = "value")

# Step 4: Plot
p <- ggplot(df_long, aes(x = iteration, y = value, color = group)) +
  geom_line(na.rm = TRUE) +
  geom_point() +
  labs(title = "Line Plot by Group", x = "Iteration", y = "Value") +
  theme_minimal()
p

# Drop the iteration column if present
cor_matrix <- cor(df_wide %>% select(-iteration), use = "pairwise.complete.obs")
cor_matrix
