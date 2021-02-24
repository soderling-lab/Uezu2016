#!/usr/bin/env Rscript

# title:
# description: 
# author: Tyler W Bradshaw

## ---- Input parameters

#BF_alpha <- 0.05 
FDR_alpha <- 0.05
FE_threshold <- 2


## ---- Set-up the workspace 

# project root dir
root <- "~/projects/Uezu2016"

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(geneLists) # for gene lists (pathways) and hyperTest function
})

# load functions in root/R and data in root/data
devtools::load_all(root, quiet = TRUE)

# load the data from root/data
data(ipsd_gene_map) # gene_map
data(surprise_partition) # partition

# Load the gene lists from twesleyb/geneLists

# geneLists TODO: rm uniprot
# geneLists TODO: enforce consistent names going forward: authorYEARword,
#   e.g. Bradshaw2021modules
# geneLists TODO: geneLists should always be associated with a publication or
# reference to additional metadata

# geneLists()

datasets <- c(# protein complexes:
	      "corum", 
	      # protein-subcellular loc annotations
	      "uniprotSubcell",
	      # ASD genes
	      "sfariGene",
	      "sfariAnimal",
	      "derubeis2014ASD", 
	      "arrabroadwave3ASD", 
	      "sanders2015ASD",
	      # other disorders
	      "alsgene",
	      "alzgene", 
	      "msgene",
	      "pdgene", # parkinsons
	      ## compiled lists, multiple disorders:
	      "g2pDBD",  # autism, epilepsy, intelllectual disability
	      "fromer2014dbd",  # combined from several studies
	      "fromer2014schz", # schz
	      "geisingerDBD", # combined ID, ADHD, BP, ect
	      ## synapse parts lists
	      "synsysnet",
	      "synGO",
	      ## epilepsy gene lists
	      "wang2017Epilepsy")

# load all requested datasets
data(list=datasets,package="geneLists")

# combine into a single list of lists
list_o_lists <- lapply(datasets, function(x) eval(parse(text=x)))
names(list_o_lists) <- datasets
gene_lists <- unlist(list_o_lists, recursive=F)

# there are a number of pathways that only contain 1 protein, 
# we wont consider these
drop <- names(which(sapply(gene_lists,length)==1))
filt_list <- gene_lists[!(names(gene_lists) %in% drop)]


## ---- Do work 

data(ipsd_bioid)

# all genes identified in BioID experiment as background
all_entrez <- unique(ipsd_bioid$Entrez)

# map partition uniprot to entrez
entrez_part <- geneLists::mapIDs(names(partition),"uniprot","entrez",gene_map)
module_list <- split(entrez_part,partition)

# rm modules with less than 2 proteins (i.e. a complex contains 2 or more proteins)
modules <- module_list[-which(sapply(module_list,length)==1)]


## ---- pre-filter pathways

# restrict to pathways with at least 1 gene in the dataset
all_paths <- names(which(sapply(filt_list, function(x) {
				   any(unlist(x) %in% all_entrez) })))


# Loop to perform hypergeometric test for each pathway and each module
results_list <- list()
pbar <- txtProgressBar(max = length(all_paths), style = 3)
for (path in all_paths) {
  # get pathway specific genes
  path_genes <- gene_lists[[path]]
  # background is all prots IDentifified
  background <- all_entrez
  # check how many pathway genes in dataset
  N <- sum(path_genes %in% background)
  # loop to perform hypergeometric test for enrichment for each module
  res_list <- list()
  for (m in seq(modules)) {
    res <- hyperTest(path_genes, modules[[m]], background)
    res[["Pathway"]] <- path
    res[["Pathway Size"]] <- length(path_genes)
    res[["Module"]] <- m
    res[["Module Size"]] <- length(modules[[m]])
    res[["n Pathway Genes"]] <- sum(modules[[m]] %in% path_genes)
    res_list[[m]] <- res
  }
  names(res_list) <- names(modules)
  # collect results in a data.table
  hyper_dt <- res_list %>% bind_rows(.id="Module")
  # adjust p-values
  hyper_dt$FDR <- p.adjust(hyper_dt$"P-value", method = "BH")
  hyper_dt$Padjust <- p.adjust(hyper_dt$"P-value", method = "bonferroni")
  # sort by fold enrichment
  hyper_dt <- hyper_dt %>% arrange(desc(`Fold enrichment`))
  # return the results
  results_list[[path]] <- hyper_dt
  setTxtProgressBar(pbar, value = match(path, all_paths))
}
close(pbar)


## ---- collect the results in a single data.table

# collect results
df <- dplyr::bind_rows(results_list) %>%
	select(Pathway, Module, `Module Size`, `Pathway Size`,
	       `n Pathway Genes`, `P-value`, FDR, Padjust, `Fold enrichment`)
df <- df %>% tidyr::separate(Pathway,into=c("Database","Pathway"),sep="\\.",extra="merge")

# only sig + enriched results:
sig_df <- df %>% filter(Padjust < FDR_alpha) %>% 
	filter(`Fold enrichment` > FE_threshold) 


## ---- save results

# save as excel
write_excel(list(GSEA=sig_df),file.path(root,"tables","iPSD_Module-GSEA.xlsx"))
message("saved :", myfile)
