#!/usr/bin/env Rscript

# title: 
# author: twab
# description: clean-up results from from running leidenalg executable

## project root dir
root <- "~/projects/Uezu2016"

## input in root/rdata
#input_part <- "surprise_partition.csv"
input_part <- "modularity_partition.csv"
#input_part <- "cpm0_partition.csv"


## ---- prepare the R env

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- import partition

# read from file
myfile <- file.path(root,"rdata",input_part)
stopifnot(file.exists(myfile))

df <- data.table::fread(myfile,drop=1) # drop first column

# munge
part <- unlist(df)

# add 1 bc python is 0-based
partition <- part + 1


## status
message("Total number of modules: ", length(unique(part))-1)


## ---- save as rda

# output filename is [input_part].rda
output_part <- paste0(tools::file_path_sans_ext(input_part),".rda")

# save partition
myfile <- file.path(root,"data", output_part)
save(partition, file = myfile, version = 2)
message("saved: ", myfile)
