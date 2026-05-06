# workflow/scripts/mit_perc_plot.R

args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]
input_dir   <- args[2]

# -------- LIBRARIES --------
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# -------- CONFIG --------
idx_dir <- input_dir

keep_chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

mt_chr <- "chrM"

# -------- CHECK INPUT --------
if (!dir.exists(idx_dir)) {
  stop("Input directory does not exist: ", idx_dir)
}

# -------- FIND FILES --------
idx_files <- list.files(
  idx_dir,
  pattern = ".idxstats.txt$",
  full.names = TRUE
)

if (length(idx_files) == 0) {
  stop("No idxstats files found in directory: ", idx_dir)
}

print(paste("Found files:", length(idx_files)))

# -------- FUNCTION --------
extract_mt_stats <- function(filepath) {
  
  df <- read_tsv(
    filepath,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  colnames(df) <- c("chrom", "length", "mapped", "unmapped")
  
  df$chrom <- str_trim(df$chrom)
  
  # TOTAL reads (only canonical chromosomes)
  df_main <- df %>% filter(chrom %in% keep_chr)
  
  total_reads <- sum(df_main$mapped, na.rm = TRUE)
  
  # MITO reads (do NOT filter them out beforehand)
  mt_reads <- sum(df$mapped[df$chrom %in% mt_chr], na.rm = TRUE)
  
  mt_fraction <- if (total_reads > 0) {
    (mt_reads / total_reads) * 100
  } else {
    NA_real_
  }
  
  tibble(
    sample = str_replace(basename(filepath), ".idxstats.txt$", ""),
    total_reads = total_reads,
    mt_reads = mt_reads,
    mt_fraction = mt_fraction
  )
}

# -------- PROCESS FILES --------
mt_stats <- bind_rows(lapply(idx_files, extract_mt_stats))

# -------- DEBUG SAFETY --------
print(mt_stats)

if (nrow(mt_stats) == 0) {
  stop("mt_stats is empty. Check input idxstats format.")
}

if (all(is.na(mt_stats$mt_fraction))) {
  stop("All mt_fraction values are NA. Check chromosome names (chrM/MT/M).")
}

# -------- WRITE CSV --------
csv_out <- file.path(idx_dir, "mitochondrial_fraction_summary.csv")
write_csv(mt_stats, csv_out)

# -------- PLOT --------
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

p <- ggplot(mt_stats, aes(
  x = mt_fraction,
  y = reorder(sample, mt_fraction)
)) +
  geom_col(fill = "steelblue") +
  theme_bw(base_size = 14) +
  labs(
    title = "Mitochondrial Reads Percentage per Sample",
    x = "Mitochondrial Reads (%)",
    y = "Sample"
  )

ggsave(output_file, plot = p, width = 10, height = 7)

print(p)