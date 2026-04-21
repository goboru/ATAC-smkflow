# workflow/scripts/plot_frip.R
args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]
input_dir   <- args[2]  # Now taking a directory path

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# List all FRiP files in the provided directory
input_files <- list.files(path = input_dir, pattern = "_frip.txt$", full.names = TRUE)

if (length(input_files) == 0) {
  stop(paste("No FRiP files found in directory:", input_dir))
}

# Read all files and combine
data_list <- lapply(input_files, function(f) {
  sample_name <- gsub("_frip.txt", "", basename(f))
  score <- as.numeric(readLines(f))
  data.frame(Sample = sample_name, FRiP = score)
})

df <- bind_rows(data_list)

# Create the plot
# Create the plot
p <- ggplot(df, aes(x = reorder(Sample, FRiP), y = FRiP, fill = FRiP)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  # Flip the plot to make long names readable on the Y-axis
  coord_flip() + 
  scale_fill_gradient(low = "firebrick3", high = "forestgreen") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "blue", alpha = 0.5) +
  theme_bw() +
  labs(title = "ATAC-seq FRiP Scores",
       subtitle = "Fraction of Reads in Peaks",
       x = "Sample Name",
       y = "FRiP Score") +
  theme(
    # Make text smaller and adjust spacing so they don't overlap
    axis.text.y = element_text(size = 7), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  ylim(0, max(df$FRiP, 0.4) + 0.05)

# Save with a wider aspect ratio to accommodate the long labels
ggsave(output_file, plot = p, width = 10, height = 8)