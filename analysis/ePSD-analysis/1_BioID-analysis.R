#!/usr/bin/env Rscript

# title: iPSD iBioID Proteomics Analysis
# description: preprocessing and statistical analysis of Dlg4 (PSD-95)
#  iBioID proteomics originally performed by AU in 2016, see Uezu2016
# author: Tyler W A Bradshaw


## ---- inputs

# project root, everything is relative to this path
root = "~/projects/Uezu2016"

# FDR and enrichment threshold 
FDR_alpha = 0.05 
enrichment_threshold = log2(4.0)

# input data in root/inst/extdata
# peptide-level data from Duke's proteomics core
input_excel = "ePSD_BioID-Supplementary-data.xlsx" 

# NOTE: Table2 has been modified to include additional sample meta_data
# NOTE: for some reason Table3 lacks protein identifier column, only contains
# 'Protein Description' column -> we map this to UniProt Accession

# raw data should be in inst/extdata
stopifnot(file.exists(file.path(root,"inst","extdata", input_excel)))


## ---- imports

suppressPackageStartupMessages({
  library(dplyr) 
  library(readxl)
  library(openxlsx)
  library(data.table) 
})

# github/twesleyb
library(getPPIs) # for mapping gene identifiers
library(tidyProt) # for statistical testing fun lmTestContrast
library(geneLists) # for a list of mito proteins

# load fun in R/ and R data in data/
# library(iPSD)
devtools::load_all(root)

# load uniprot to gene identifier mapping data
data(list="uniprot_map",package="geneLists")

# load mitochondrial protein list from twesleyb/geneLists
data(list="mitocarta2", package = "geneLists")


## ---- Load the raw proteomics data 

# raw excel data is in root/inst/extdata
data_list <- list()
myfile <- file.path(root,"inst","extdata", input_excel)

# loop to read excel sheets
excel_sheets <- readxl::excel_sheets(myfile)
for (sheet in excel_sheets) {
	data_list[[sheet]] <- readxl::read_xlsx(myfile, sheet = sheet)
}

# MS info
ms_info <- data_list[[1]] 

# sample meta data
# NOTE: user must edit this table
meta_data <- data_list[[2]] %>% 
	dplyr::select(Sample,Condition,Subject,BioReplicate,Mixture)

# raw, peptide-level data
raw_pep <- data_list[[3]] 

# load sheet 4 (gene map)
gene_map <- data_list[[4]] %>% 
	tidyr::separate('Protein Names', into=c("ID","Accession"), sep=", ") %>% 
	mutate(Accession=gsub("SP\\:","",Accession))

# finish mapping uniprot to entrez and gene symbols
gene_map$uniprot <- gsub("\\-*","", gene_map$Accession)

# method 1: using uniprot_map
entrez <- geneLists::mapIDs(gene_map$uniprot,
			    from="Accession",to="Entrez",
			    gene_map=uniprot_map)

# replace empty string with NA
entrez[entrez==""] <- NA

warning("Unable to map ",sum(is.na(entrez))," UniProt to Entrez identifiers.") # 215

# there are many prots with missing entrez and symbol bc its uniprot accession is
# outdated!

# ignore wanrings (uniprot not mapped to entrez)
gene_map <- gene_map %>% 
	mutate(entrez) %>% 
	mutate(symbol=getIDs(entrez,from="entrez",to="symbol",species="mouse"))

# NOTE: several hundred uniprot IDs including psd95's are not mapped
# To fix these edge cases (resulting from outdated uniprot accession),
# we extract gene Names from 'Protein Description' column.

# then use these gene names to get entrez IDs associated with each gene

# munge to extract gene symbol from GN= in description
d <- gene_map$"Protein Description"
namen <- gsub("GN=","", sapply(strsplit(regmatches(d, regexpr("GN=.*",d)),"\\ "),"[",1))
idx <- grepl("GN=.*",d)
names(namen) <- d[idx]
gene_map$GN <- namen[gene_map$"Protein Description"]

# map these symbols to entrez
# ignore warnings
gene_map <- gene_map %>% mutate(entrez2 = getIDs(GN,"symbol","entrez","mouse"))

# fix missing entrez and symbols
idx <- is.na(gene_map$entrez) & !is.na(gene_map$entrez2)
gene_map$entrez[idx] <- gene_map$entrez2[idx]
gene_map$symbol[idx] <- gene_map$GN[idx]

# there are a small number of genes with missing gene identifier info
# several of these proteins appear to be contaminants
idx <- is.na(gene_map$entrez)

# final  attempt to map using uniprot IDs
gene_map$entrez[idx] <- mapIDs(gene_map$ID[idx],"ID","Entrez",uniprot_map)

# drop any remaining prot with missing entrez info
to_drop <- gene_map$uniprot[is.na(gene_map$entrez)]

# method 2: using getIDs (org.mm.eg.db)
#entrez <- geneLists::getIDs(uniprot,from="uniprot",to="entrez",species="mouse")
#sum(is.na(entrez))

# method3 3: using mgi_map
#data(mgi_map, package="geneLists")
#entrez <- geneLists::mapIDs(uniprot,from="UniProt",to="Entrez",gene_map=mgi_map)
#sum(is.na(entrez))


## ---- annot raw protein data with protein identifiers

idx <- match(raw_pep$"Protein Description", gene_map$"Protein Description")
raw_pep <- raw_pep %>% mutate(Protein = gene_map$uniprot[idx]) %>% 
	# remove small number of unmapped proteins
	filter(!(Protein %in% to_drop))

# check initial number of prots
prots <- unique(raw_pep$Protein)
knitr::kable(data.table(nProts=length(prots)))

stopifnot(!any(is.na(prots)))


## ---- tidy-up the raw data

# clean-up column names for QC
# columns with data are distinguished by "_"
# id columns are all the cols that DONT contain data:
id_cols <- colnames(raw_pep)[!grepl("_",colnames(raw_pep))]

# tidy the data, melt into long format
tidy_pep <- raw_pep %>% 
	# melt into long format
	reshape2::melt(id.var = id_cols, 
		       variable.name="Sample", 
		       value.name = "Intensity", 
		       variable.factor = FALSE) %>%
        mutate(Sample=as.character(Sample))

samples <- unique(tidy_pep$Sample)
stopifnot(all(samples %in% meta_data$Sample))

# merge with sample meta data
tidy_pep <- tidy_pep %>% left_join(meta_data, by = "Sample")


## ----- remove mitochondrial contaiminants

# collect mitocarta entrez ids
mito_entrez <- unlist(mitocarta2, use.names=FALSE)

# rm these additional mito prots
hand_anno_mito <- c("Mtres1", "Gcdh", "Clpx", "Pdk3", "Mrps36","Hscb",
	  "Aldh2","Shmt2","Ciapin1","Ssbp1","Bckdha","Dap3")
hand_anno_entrez <- geneLists::getIDs(hand_anno_mito, 
				      'symbol', 'entrez', 'mouse')

# to be removed if in mito
mito <- c(mito_entrez,hand_anno_entrez) 
mito_prot <- mapIDs(mito,"Entrez","Accession",gene_map=uniprot_map)

# total number of mito prots to be removed
nMito <- sum(mito_prot %in% tidy_pep$Protein)
message("Mitochondrial proteins removed: ", nMito)

# rm mito prots
tidy_pep <- tidy_pep %>% dplyr::filter(!(Protein %in% mito_prot))
tidy_pep$Sample <- as.character(tidy_pep$Sample)
tidy_pep$Condition <- as.character(tidy_pep$Condition)


## ---- sample loading normalization 

# sl normalization
sl_pep <- normSL(tidy_pep, groupBy="Sample")

# Check, column sums (Run total intensity) should now be equal:
sl_pep %>% 
	group_by(Sample) %>% 
	summarize("Total Intensity"=sum(Intensity,na.rm=TRUE),.groups="drop") %>%
	knitr::kable()


## ---- protein summarization 

# protein summarization by sum
raw_prot <- sl_pep %>% group_by(Protein, Sample) %>% 
	summarize(Peptides = sum(!is.na(Intensity)), 
		  Intensity = sum(Intensity,na.rm=TRUE),.groups="drop") %>%
        left_join(meta_data,by="Sample")


## ---- protein level filtering 

# Remove one hit wonders -- proteins identified by a single peptide
one_hit_wonders <- unique(raw_prot$Protein[raw_prot$Peptides == 1])
message("Number of one-hit-wonders: ", length(one_hit_wonders))


## ---- sample pool normalization 

# perform normalization to QC samples
spn_prot <- normSP2(raw_prot, pool="QC")


## ---- statistical testing for positive control, e.g. psd95

# get psd95's uniprot
idx <- grepl("Dlg4", gene_map$"Protein Description")
psd95 <- gene_map$uniprot[idx]
stopifnot(psd95 %in% spn_prot$Protein)

# simple linear model:
fx <- log2(Intensity) ~ 0 + Condition

# * NOTE: by setting the intercept to 0 we explicitly estimate 
# all levels of Condition

# * NOTE: you may alternatively set QC as reference level:
# fx <- log2(Intensity) ~ Condition | levels(Condition)[1] == "SPQC"

# fit the model 
# remove QC data prior to testing (var = 0)
df <- spn_prot %>% subset(Protein == psd95) %>% filter(Condition != "SPQC")
fm <- lm(fx, data=df)

# create a contrast:
LT <- tidyProt::getContrast(fm,"DLG4-BioID","Control")

lmTestContrast(fm, LT) %>% mutate(Pvalue=formatC(Pvalue)) %>% knitr::kable()


## ---- need complete cases!

# rm 0
spn_prot$Intensity[spn_prot$Intensity == 0] <- NA

# remove proteins with any remaining mising values
dm <- spn_prot %>% dcast(Protein ~ Sample, value.var = "Intensity")
idx <- apply(dm,1,function(x) any(is.na(x)))
to_rm <- dm$Protein[idx]

message("drop proteins with missing values: ", length(to_rm))
norm_prot <- spn_prot %>% filter(!(Protein %in% to_rm))

# final number of testable proteins
knitr::kable(data.table(nProts=formatC(length(unique(norm_prot$Protein)))))


## ---- loop to perform test for all proteins

# empty lists for models, sigma2, and degrees of freedom for each fit
fm_list <- list()
s2_list <- list()
df_list <- list()

# fit the models
proteins <- unique(norm_prot$Protein)
for (prot in proteins) {
	df <- norm_prot %>% 
		subset(Protein == prot) %>% 
		filter(Condition != "SPQC")
	fm <- lm(fx, data = df)
	fm_list[[prot]] <- fm
	s2_list[[prot]] <- sigma(fm)^2 # calculate sigma2 for moderation
	df_list[[prot]] <- fm$df.residual # get df for moderation
}

# perform t-statistic moderation using limma
eb_var <- limma::squeezeVar(unlist(s2_list), unlist(df_list))
df_prior <- eb_var$df.prior
s2_prior <- eb_var$s2.prior
if (is.null(s2_prior)) { s2_prior <- 0 }

# s2 prior == 0, interpretation??
knitr::kable(data.table(DF=df_prior,S2=s2_prior))

# examine gphn result with moderation
df <- norm_prot %>% subset(Protein == psd95) %>% filter(Condition != "SPQC")
fm <- lm(fx, data=df)

# result with moderation
lmTestContrast(fm, LT, s2_prior, df_prior) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% knitr::kable()

# loop to moderate comparisons for all proteins
result_list <- list()
for (prot in proteins) {
	fm <- fm_list[[prot]]
	result_list[[prot]] <- lmTestContrast(fm,LT,s2_prior,df_prior)
}

# collect results
results_df <- dplyr::bind_rows(result_list,.id="Protein")

# calculate FDR for each contrast
results_df <- results_df %>% group_by(Contrast) %>% 
	mutate(FDR = p.adjust(Pvalue, method="BH"))

# anno with sig and up
results_df <- results_df %>% 
	mutate(up = log2FC > enrichment_threshold) %>% 
	mutate(sig = FDR < FDR_alpha) %>%
	mutate(candidate = up & sig)

# annotate with gene Symbols and Entrez IDs
results_df <- results_df %>%
       	mutate(Entrez = mapIDs(Protein,"uniprot","entrez",gene_map)) %>%
	mutate(Symbol = getIDs(Entrez,"entrez","symbol","mouse"))


# add GN= annotation
#results_df <- results_df %>%
#	mutate(GeneName=mapIDs(Protein,"uniprot","GN",gene_map))

# final sort
results_df <- results_df %>% 
	arrange(desc(sig), desc(up), Pvalue, desc(log2FC)) %>%
	dplyr::select(-sig,-up)

# summary for each contrast
results_df %>% 
	group_by(Contrast) %>% 
	summarize(nSig=sum(candidate),.groups="drop") %>% 
	knitr::kable()

# sort cols
results_df <- results_df %>% dplyr::select(Protein,Entrez,Symbol,Contrast,log2FC,percentControl,SE,Tstatistic,Pvalue,DF,FDR,candidate)

# check psd95
results_df %>% filter(Protein == psd95) %>% knitr::kable()


## ---- save results

# annot normalized protein data with symbols and entrez
spn_prot <- spn_prot %>% 
	mutate(Symbol = mapIDs(Protein,"uniprot","symbol",gene_map)) %>%
	mutate(Entrez = mapIDs(Protein,"uniprot","entrez",gene_map))

# final munge...

results_list <- results_df %>% group_by(Contrast) %>% group_split()

namen <- gsub("Condition|Control|-","",
	      sapply(results_list, function(x) unique(x$Contrast)))
names(results_list) <- namen

class(results_list) <- "list"

results_list[["Raw Data"]] <- as.data.frame(raw_pep)
results_list[["Norm Data"]] <- spn_prot

# save tidy, normalized protein data as rda
myfile <- file.path(root,"data","epsd_bioid.rda") 
epsd_bioid <- spn_prot
save(epsd_bioid,file=myfile,version=2)
message("saved: ", myfile)

# save statistical results as excel
myfile <- file.path(root,"tables","ePSD_iBioID_Results.xlsx")
write_excel(results_list, myfile)
message("saved: ", myfile)

# save statistical results 
myfile <- file.path(root,"data","epsd_results.rda")
epsd_results <- results_df
save(epsd_results,file=myfile,version=2)
message("saved: ", myfile)

# save gene map
myfile <- file.path(root,"data","epsd_gene_map.rda")
save(gene_map,file=myfile,version=2)
message("saved: ", myfile)
