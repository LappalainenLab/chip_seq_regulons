#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(stringr)
library(argparser)
library(here)

here::i_am("README.md")

# Set working directory
setwd(here())

# Parse command-line arguments
# Create a parser
p <- arg_parser("Get input file")

# Add command line arguments
p <- add_argument(p, "file_name", help="data file name", type="character")
argv <- parse_args(p)

# Read TF target mapping data
data = fread(paste0(argv$file_name), nThread=5)

lookup <- c(distance_S2 = "distance.x", distance_M2 = "distance.y")
rename(data, all_of(lookup)) -> data

data %>% select(-c(color_start.x, color_end.x, color.x, color_start.y, color_end.y, color.y)) %>% distinct() -> data

fwrite(data, paste0(str_split_1(argv$file_name, pattern="[.]")[1], "_cleaned.tsv"), col.names=T, row.names=F, quote=F, sep="\t")
