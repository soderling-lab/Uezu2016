#!/usr/bin/env Rscript

# title: iPSD Proteomics
# author: twab
# description: compiling ePSD PPIs

## ---- inputs

root <- "~/projects/iPSD"


## ---- imports

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(data.table)
})

library(getPPIs)
data(musInteractome,package="getPPIs")

# library(iPSD)
devtools::load_all(root)

data(epsd)
data(epsd_gene_map)



## ---- let's test a hypothesis

library(geneLists)

geneLists()

data(sfariGene,package="geneLists")
data(sfariAnimal,package="geneLists")

# this is a list of genes
str(sfariAnimal)

# entrez IDs are the most stable gene identifier
sfari <- c(sfariGene[["ASD"]], sfariAnimal[["ASD"]])


data.table(`SFARI ASD Genes` = length(sfari)) %>% knitr::kable()


sum(epsd %in% sfari)
length(epsd)

# test for enrichment
data(epsd_bioid)

all_genes <- unique(epsd_bioid$Entrez)

# approximate two fold enrichment of SFARI genes in ePSD
geneLists::hyperTest(sfari,epsd,background=all_genes) %>% t() %>% knitr::kable()

## ---- main

# use geneLists::getIDs to map entrez to gene symbols
symbols <- toupper(geneLists::getIDs(entrez,"entrez","symbol","mouse"))

noa_df <- data.table(entrez=epsd,
		     symbol=symbols) %>% 
            mutate(SFARI= entrez %in% sfari)

# write to file
fwrite(noa_df,file="noa_df.csv")

## ---- collect ppis (edges)

entrez <- epsd

os_keep <- c(9606, 10116, 10090)

ppi_df <- musInteractome %>% 
	filter(Interactor_A_Taxonomy %in% os_keep) %>%
	filter(osEntrezA %in% entrez) %>% 
	filter(osEntrezB %in% entrez) %>% 
	dplyr::select(osEntrezA,osEntrezB,Publications)


idx <- match(ppi_df$osEntrezA,gene_map$entrez)
idy <- match(ppi_df$osEntrezB,gene_map$entrez)
ppi_df <- ppi_df %>% mutate(protA = gene_map$uniprot[idx],
		  protB = gene_map$uniprot[idy])


## ---- convert ppi_df to adjm

g <- ppi_df %>% 
	dplyr::select(protA,protB) %>% 
	graph_from_data_frame(directed=FALSE) %>% 
	simplify()

adjm <- as.matrix(as_adjacency_matrix(g))

## TODO: examine topology!!!


## ---- save data

# save ppi adjm for clustering
myfile <- file.path(root,"rdata","ppi_adjm.csv")
fwrite(adjm %>% as.data.table(keep.rownames="protein"), myfile)

# save node attributes
myfile <- file.path(root,"rdata","epsd_noa.csv")
fwrite(noa_df,myfile)
message("saved: ", myfile)

# save edge data.frame
myfile <- file.path(root,"rdata","epsd_ppi.csv")
fwrite(ppi_df,myfile)
message("saved: ", myfile)
