#!/usr/bin/env Rscript

## ---- inputs

root <- "~/projects/Uezu2016"


## ---- imports

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  # for working with graphs:
  library(igraph)
})

# twesleyb/geneLists
library(geneLists)

# load mouse PPI dataset
data(musInteractome,package="getPPIs")

# load project's data
devtools::load_all(root)

data(insyn1)
data(ipsd_gene_map) # gene_map
data(ipsd_sig_prots) # sig_prots
data(ipsd_colors) # ipsd_colors == colors for sig_prots


## ---- load other partitions

parts <- c("surprise", 
	   "significance", 
	   "cpm0", 
	   "cpm1", 
	   "rber0", 
	   "rber1", 
	   "rbconfiguration0", 
	   "rbconfiguration1", 
	   "modularity")

# collect filepaths to partition files in root/rdata
parts_files <- sapply(parts, function(x) {
			      list.files(file.path(root,"rdata"),
					 pattern=x, 
					 full.names=T)
	   })

# read from file
part_list <- lapply(lapply(parts_files, fread, drop=1),unlist)

# set small modules to 0
part_list <- lapply(part_list, function(p) { 
			    p[p %in% as.numeric(names(which(table(p)<3)))] <- 0
			    return(p) 
	   })

# examine total number of modules for each method
k <- sapply(part_list,function(x) length(unique(x)))
k %>% t() %>% knitr::kable()


## ---- main

# create node attribute file for all sig_prots
prots <- unlist(sig_prots)
idx <- match(prots,gene_map$uniprot)
noa_df <- data.table(protein=prots, 
		     gene=gene_map$symbol[idx], 
		     entrez=gene_map$entrez[idx]) %>% 
  mutate(SYMBOL=toupper(gene))

# simplify sig_prot names
names(sig_prots) <- gsub("Condition|BioID|Control|-","", names(sig_prots))

# collect node attributes
noa_df <- noa_df %>% 
	mutate(cbBioID = protein %in% sig_prots[["Collybistin"]]) %>%
	mutate(isBioID = protein %in% sig_prots[["InSyn1"]]) %>%
	mutate(gpBioID = protein %in% sig_prots[["Gephyrin"]]) %>%
	mutate(overlap = protein %in% sig_prots[["overlap"]])


## ---- map sig_prots to color annotation

# ipsd_colors
#|Gphn    |Arhgef9 |InSyn1  |Gphn+Arhgef9 |Gphn+Insyn1 |Arhgef9+Insyn1 |All    |
#|:-------|:-------|:-------|:--------- --|:-----------|:--------------|:------|
#|#4CFF59 |#FEFF5A |#4C3EFA |#BDFF55      |#44BF81     |#7F7CBB        |#85D979|

prots <- noa_df$protein
colors <- setNames(rep(col2hex("gray"),length(prots)), nm=prots)

idx_A <- names(colors) %in% sig_prots[["Gephyrin"]]
idx_B <- names(colors) %in% sig_prots[["Collybistin"]]
idx_C <- names(colors) %in% sig_prots[["InSyn1"]]

colors[idx_A] <- "#4CFF59"
colors[idx_B] <- "#FEFF5A"
colors[idx_C] <- "#4C3EFA" 
#
idx_AB <- idx_A & idx_B
colors[idx_AB] <- "#BDFF55"
#
idx_AC <- idx_A & idx_C
colors[idx_AC] <- "#44BF81"
#
idx_BC <- idx_B & idx_C
colors[idx_BC] <- "#7F7CBB"
#
idx_ABC <- idx_A & idx_B & idx_C
colors[idx_ABC] <- "#85D979"
#
# annotate with colors
noa_df$color <- colors[noa_df$protein]


## ---- collect ppis (edges)

entrez <- noa_df$entrez

# mouse human and rat
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

message("Compiled ", dim(ppi_df)[1], " PPIs.")


## ---- add select interactions tested by Uezu2016

iqsec3 <- mapIDs("Iqsec3","symbol","uniprot",gene_map)
arhgap32 <- mapIDs("Arhgap32","symbol","uniprot",gene_map)

edges <- list(c(iqsec3,gphn), # IQSEC3 <-> GPHN
	      c(insyn1,gphn), # INSYN1 <-> GPHN
	      c(insyn2,gphn), # INSYN2 <-> GPHN
	      c(arhgap32,gphn)) # ARHGAP32 <-> GPHN

my_edges <- lapply(edges,function(x) 
       data.table(osEntrezA = mapIDs(x[1],"uniprot","entrez",gene_map),
		  osEntrezB = mapIDs(x[2],"uniprot","entrez",gene_map),
		  Publications="pubmed:27609886",protA=x[1],protB=x[2]))

# add to edge data
ppi_df <- rbind(ppi_df, bind_rows(my_edges))


## ---- convert ppi_df to adjm

g <- ppi_df %>% 
	dplyr::select(protA,protB) %>% 
	graph_from_data_frame(directed=FALSE) %>% 
	simplify()

# convert to adjm with igraph's as_adjacency_matrix
adjm <- as.matrix(as_adjacency_matrix(g))


## ---- examine toplogy

# connectivity histogram
#hist(apply(adjm,1,sum))

# connectivity plot
plot <- ggplotScaleFreeFit(apply(adjm,1,sum))

# save plot
myfile <- file.path(root,"figs","scale-free-fit.pdf")
ggsave(plot,file=myfile,width=3.0,height=3.0)
message("saved: ", myfile)


## insure noa_df has membership info

# surprise partition
noa_df$membership <- as.numeric(partition[noa_df$protein])
noa_df$membership[is.na(noa_df$membership)] <- 0

names(part_list)
noa_df$cpm1 <- part_list[["cpm1"]][noa_df$protein]
noa_df$cpm1[is.na(noa_df$cpm1)] <- 0

data(surprise_colors)
noa_df$surprise <- module_colors[paste0("M",noa_df$membership)]
noa_df$surprise[is.na(noa_df$surprise)] <- col2hex("gray")

data(modularity_colors)
data(modularity_partition)
m <- partition[noa_df$protein]
m[is.na(m)] <- 0
noa_df$modularity <- module_colors[paste0("M",m)]
noa_df$modularity[is.na(noa_df$modularity)] <- col2hex("gray")

data(cpm1_colors)
data(cpm1_partition)
m <- partition[noa_df$protein]
m[is.na(m)] <- 0
noa_df$cpm1 <- module_colors[paste0("M",m)]
noa_df$cpm1[is.na(noa_df$cpm1)] <- col2hex("gray")


## ---- save data

# save ppi adjm for clustering
myfile <- file.path(root,"rdata","ppi_adjm.csv")
fwrite(adjm %>% as.data.table(keep.rownames="entrez"), myfile)

# save node attributes
myfile <- file.path(root,"rdata","noa_df.csv")
fwrite(noa_df,myfile)
message("saved: ", myfile)

# save edges
myfile <- file.path(root,"rdata","ppi_df.csv")
fwrite(ppi_df,myfile)
message("saved: ", myfile)

# import these files into Cytoscape
