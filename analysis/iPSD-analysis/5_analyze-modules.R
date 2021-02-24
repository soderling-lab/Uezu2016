#!/usr/bin/env Rscript

# project root
root <- "~/projects/iPSD"

# library(iPSD)
devtools::load_all(root)

# load the partition
data(ipsd_gene_map)
data(ipsd_partition)

# imports
suppressPackageStartupMessages({ 
	library(dplyr)
	library(data.table)
})

uniprot <- names(partition)
idx <- match(uniprot,gene_map$uniprot)

#symbols = geneLists::getIDs(entrez,from="entrez",to="symbol","mouse")

df <- data.table(symbol=gene_map$symbol[idx],
		 protein=gene_map$uniprot[idx],
		 membership=partition)


data(ipsd_sig_prots)

# > names(sig_prots)
# [1] "Collybistin" "Gephyrin"    "InSyn1"      "overlap"    
names(sig_prots) <- gsub("Condition|BioID|Control|-","", names(sig_prots))

df <- df %>% mutate(Collybistin = protein %in% sig_prots[["Collybistin"]],
	            Insyn1 = protein %in% sig_prots[["InSyn1"]],
	            Gephyrin = protein %in% sig_prots[["Gephryin"]],
	            Common = protein %in% sig_prots[["overlap"]])

# add entrez annot
df <- df %>% mutate(entrez = geneLists::mapIDs(protein,"uniprot","entrez",gene_map))

# m12 broken down by BioID bait
df %>% filter(membership==12) %>% summarize(sum(Collybistin),sum(Insyn1),sum(Gephyrin),sum(Common))

## compare with corum
#table(partition)
#
#
#data(corum,package="geneLists")
#
## list of mouse entrez genes
## names(corum)
## x = corum[[1]]
## geneLists::getIDs(x,"entrez","symbol","mouse")
#
#data(ipsd_results)
#ipsd <- unique(ipsd_results$Entrez[ipsd_results$candidate])
#
#library(geneLists)
#
#module_list <- split(mapIDs(names(partition),"uniprot","entrez",gene_map),partition)
#gene_list <- corum
#
#data(ipsd_bioid)
#
#all_genes <- unique(ipsd_bioid$Entrez)
#
#idx <- grepl("Dystrophin",names(corum))
#sum(idx)
#names(corum)[idx]
#
#dystrophin_complex <- unique(unlist(corum[idx],use.names=F))
#
#
#m12 <- mapIDs(names(partition)[partition==12],"uniprot","entrez",gene_map)
#
#sum(m12 %in% dystrophin_complex)
#
#hyperTest(dystrophin_complex,m12,background=all_genes)
#
#df %>% filter(symbol=="Dmd")
#df %>% filter(membership==12) 
