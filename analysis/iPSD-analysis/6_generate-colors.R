#!/usr/bin/env Rscript

# title: 
# description: generate module colors
# authors: Tyler W Bradshaw

## ---- INPUT

# input data in root/data
root = "~/projects/Uezu2016"

part_file = "surprise_partition"

# not clustered == "gray"
NC_color = "#BEBEBE" 

## ---- OUTPUT 
# * color assignments for every module in graph partition


## ---- FUNCTIONS 

str_to_vec <- function(response) {
  # parse the python dictionary returned as a string from 
  # system(random_color.py)
	vec <- gsub("'","",gsub("\\[|\\]","",
				trimws(unlist(strsplit(response,",")))))
	return(vec)
}


## ---- Prepare the workspace 

# Global Imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Local Imports 
devtools::load_all(root, quiet=TRUE)

# Load TMT data and partition
data(list=part_file)
message("loaded: ", part_file)

# rm small modules
#min_size = 3
#idx <- partition %in% as.numeric(names(which(table(partition)<min_size)))
#partition[idx] <- 0


## ---- Generate Module colors 

prots <- names(partition)

modules <- split(prots, partition)
names(modules) <- paste0("M",names(modules))

n_colors <- length(modules) 
message("generating: ", n_colors, " random colors")

# Path to python script which is a simple script that uses the python 
# port of randomcolors to generate random colors.
script <- file.path(root,"Py","random_color.py")

stopifnot(file.exists(script))

# how can we conveniently package this python code for R?

# Generate n random colors
cmd <- paste(script,"--count", n_colors, "--luminosity", "bright")
response <- system(cmd, intern = TRUE)

# NOTE: to install randomcolor, use conda:
# conda install -c conda-forge randomcolor

#  Parse the response
colors <- toupper(str_to_vec(response))

# Module color assignments
# Initialize a vector for the module colors
module_colors <- setNames(rep(NA, length(modules)), nm=names(modules))

# Insure that M0 is gray 
module_colors["M0"] <- NC_color

# The reamining colors are random
idx <- is.na(module_colors)
module_colors[idx] <- sample(colors,sum(idx))


## ---- Save the data 

# Save updated module colors
namen <- paste0(gsub("partition","colors",part_file),".rda")
myfile <- file.path(root,"data", namen)
save(module_colors,file=myfile,version=2)
message("saved: ", myfile)
