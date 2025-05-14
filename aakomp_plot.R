#!/usr/bin/env Rscript

# aaKomp CDF plot script
# Written by Johnathan Wong

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript aakomp_plot.R <input_file> <target_value> <score_file>\n")
  quit(status = 1)
}

in_file <- args[1]
target_line <- as.numeric(args[2])
score_file <- args[3]
base <- basename(in_file)

# Read the score value from the score file
score <- tryCatch({
  as.numeric(readLines(score_file, warn = FALSE)[1])
}, error = function(e) {
  cat(sprintf("Could not read score from %s\n", score_file))
  quit(status = 1)
})

if (is.na(score)) {
  cat(sprintf("Invalid numeric score in %s\n", score_file))
  quit(status = 1)
}

# Read input file (skipping header)
data <- suppressWarnings(read_delim(
  in_file,
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE,
  skip = 1, show_col_types = FALSE
))

if (ncol(data) < 9) {
  cat(sprintf("File %s does not have enough columns\n", base))
  quit(status = 1)
}

colnames(data)[6] <- "ratio"
colnames(data)[9] <- "id"

vals <- data %>%
  filter(!is.na(ratio)) %>%
  group_by(id) %>%
  summarise(max_ratio = max(ratio), .groups = "drop") %>%
  arrange(max_ratio) %>%
  mutate(cdf = row_number() / n())

if (nrow(vals) == 0) {
  cat(sprintf("No valid entries in %s\n", base))
  quit(status = 1)
}

plot <- ggplot(vals, aes(x = max_ratio, y = cdf)) +
  geom_step(direction = "hv", linewidth = 1) +
  geom_vline(xintercept = target_line, linetype = "dotted", color = "red", linewidth = 1) +
  annotate("text", x = target_line, y = 0.2, angle = 90, vjust = -0.5,
           label = "Targeted miBF rescue threshold", color = "red", size = 5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("aaKomp score = %.2f", score), size = 5) +
  labs(
    title = paste("aaKomp Cumulative Distribution Function (CDF) –", base),
    x = expression(italic(k)*"-mer completeness ratio"),
    y = "CDF"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20)
  )

# Save the plot
ggsave(filename = sprintf("%s_cdf.png", in_file), plot = plot, width = 13.5, height = 9, dpi = 300)
