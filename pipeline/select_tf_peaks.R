#!/usr/bin/env Rscript

library(data.table)
library(argparser)
library(dplyr)
library(here)

here::i_am("README.md")

# Set working directory
setwd(here())

# Create a parser
p <- arg_parser("Get input file")

# Add command line arguments
p <- add_argument(p, "tf", help="tf name", type="character")
p <- add_argument(p, "file", help="peak file path", type="character")
p <- add_argument(p, "path", help="output path", type="character")

# Parse the command line arguments
argv <- parse_args(p)

print(argv)

data = fread(argv$file, nThread=10)

data %>% filter(is_S2Mb | is_M2Kb | is_S2Kb) %>%
        filter(tf == toupper(argv$tf)) -> tss_bed

fwrite(tss_bed, paste0(argv$path, "/", argv$tf, "_temp_peak_file.bed"), quote=F, sep="\t", row.names = F, col.names = F)





