#!/usr/bin/env Rscript

# title: iPSD iBioID Proteomics Analysis
# description: preprocessing and statistical analysis of InSyn1, Gephyrin, and
# Collybistin iBioID proteomics originally performed by AU in 2016, see Uezu2016
# author: Tyler W A Bradshaw





## ---- inputs

# project root, everything is relative to this path
root = "~/projects/soderling-lab/Uezu2016"

# FDR and enrichment threshold 
FDR_alpha = 0.05 
enrichment_threshold = log2(4.0)

# input data in root/inst/extdata
# peptide-level data from Duke's proteomics core
input_excel = "iPSD_BioID-Supplementary-data.xlsx" 

# should be in inst/extdata
stopifnot(file.exists(file.path(root,"inst","extdata", input_excel)))


## ---- imports

suppressPackageStartupMessages({
  library(dplyr) 
  library(readxl)
  library(openxlsx)
  library(data.table) 
})

# R > devtools::install_github(twesleyb/*)
# * library(getPPIs) # for mapping gene identifiers
library(tidyProt) # for statistical testing fun lmTestContrast
library(geneLists) # for a list of mito proteins and gene mapping functions

# load fun in R/ and R data in data/
devtools::load_all(root)

# load uniprot to gene identifier mapping data
data(list="uniprot_map",package="geneLists")

# load mitochondrial protein list from twesleyb/geneLists
data(list="mitocarta2", package = "geneLists")


## ---- other functions

normSP <- function(tp, pool){
	# perform normalization to pooled QC samples
	# store a copy of the data
	tp <- ungroup(tp)
	tp_copy <- tp
	# group pooled together
	tp$Group <- as.numeric(grepl(paste(pool,collapse="|"),tp$Sample))
	tp_list <- tp %>% group_by(Protein,Group) %>% 
		dplyr::summarize(Mean_Intensity=mean(Intensity,na.rm=TRUE),
			  n = length(Intensity), .groups="drop") %>%
	as.data.table() %>% 
	arrange(Protein,Group) %>% 
	group_by(Protein) %>% group_split()
	# loop to calculate normalization factors
        new_list <- list()
	for (i in 1:length(tp_list)){
		x <- tp_list[[i]]
		x$NormFactor <- c(x$Mean_Intensity[2]/x$Mean_Intensity[1],1)
		x$Norm_Mean_Intensity <- x$Mean_Intensity * x$NormFactor
		new_list[[i]] <- x
	}
	tp_list <- new_list
	# collect in a df
	df <- do.call(rbind,tp_list) %>% 
		dplyr::select(Protein,Group,NormFactor)
	# merge with input data
	tp_norm <- left_join(tp,df,by=c("Protein","Group"))
	# Perform normalization step.
	tp_norm$Intensity <- tp_norm$Intensity * tp_norm$NormFactor
	tp_norm <- tp_norm %>% dplyr::select(colnames(tp_copy))
	tp_norm <- as.data.table(tp_norm)
	# return the normalized data
	return(tp_norm)
} #EOF


## ---- Load the raw proteomics data 

data_list <- list()
myfile <- file.path(root,"inst","extdata", input_excel)

# loop to read excel sheets
excel_sheets <- readxl::excel_sheets(myfile)
for (sheet in excel_sheets) {
	data_list[[sheet]] <- readxl::read_xlsx(myfile, sheet = sheet)
}

# names(data_list)
# [1] "Table 1" "Table 2" "Table 3"

# MS info - we don't need this
ms_info <- data_list[[1]] 

# sample meta data - we added column Condition!
meta_data <- data_list[[2]] %>% 
	dplyr::select(Sample,Condition)

# raw, peptide-level data - our starting data
raw_pep <- data_list[[3]] 


## ---- filter raw protein data - rm non-mouse species

# we need to: 
# [1] extract Uniprot identifiers from 'Primary Protein Name' column
# [2] map this to something more stable, entrez IDs are best

# protein ids are uniprot KB ids!
identifier_type <- "UniProtID"
protein_identifier_col <- "Primary Protein Name"

# head(sample(all_ids))
# [1] "ATPG_MOUSE"  "IF4G1_MOUSE" "SCAM3_MOUSE" "PTPA_MOUSE"  "TAB2_MOUSE" 
# [6] "AGK_MOUSE"  
all_ids <- unique(raw_pep[[protein_identifier_col]])

# Unable to map 142 input protein identifiers. 
# ^ Uniprot IDs are not great identifiers bc they are very dynamic
not_mapped <- is.na(match(all_ids,uniprot_map$"UniProtKB-ID"))
warning("Unable to map ", sum(not_mapped), " input protein identifiers.")

## make sure insyn1 is in the data
idx <- uniprot_map$GeneID == 319477
insyn1 <- uniprot_map$"UniProtKB-AC"[idx]

## make sure insyn2a is in the data
idx <- uniprot_map$GeneID == 627214
insyn2 <- uniprot_map$"UniProtKB-AC"[idx]

# map uniprot ids to accession
uniprot <- uniprot_map$"UniProtKB-AC"[match(all_ids, uniprot_map$"UniProtKB-ID")]
uniprot <- setNames(uniprot,nm=all_ids)

# get insyn1's name and map to uniprot
idx <- grepl("Uncharacterized protein C3orf67 homolog*", 
	     raw_pep$"Protein Description")
namen <- unique(raw_pep[[protein_identifier_col]][idx])
uniprot[namen] <- insyn1

# get insyn2's name and map to uniprot
idx <- grepl("Fam196a", raw_pep$"Protein Description")
namen <- unique(raw_pep[[protein_identifier_col]][idx])
uniprot[namen] <- insyn2

# map uniprot accession to entrez
entrez <- geneLists::queryMGI(uniprot)

# insure insyn2 entrez is correct
entrez <- setNames(entrez, nm=uniprot)
entrez[insyn2] <- 627214

# combine as gene_map
gene_map  <- data.table(id=all_ids, uniprot=uniprot, entrez=entrez) %>%
	filter(!is.na(entrez)) %>%
	mutate(symbol = getIDs(entrez,"entrez","symbol","mouse"))

# we can use gene_map to look-up gene identifiers

# insure that we got insyn1
stopifnot("Insyn1" %in% gene_map$symbol)
stopifnot("Insyn2a" %in% gene_map$symbol)


## ---- tidy-up the raw data

# clean-up column names for QC
# id columns are all the cols that DONT contain data:
id_cols <- colnames(raw_pep)[!grepl("-",colnames(raw_pep))]

# drop unmapped proteins (by inspection these are
# typically contaminants or edge case genes)
idx <- is.na(match(raw_pep[[protein_identifier_col]],gene_map$id))
# tidy the data, melt into long format
tidy_pep <- raw_pep %>% 
	filter(!idx) %>%
	# melt into long format
	reshape2::melt(id.var = id_cols, 
		       variable.name="Sample", 
		       value.name = "Intensity", 
		       variable.factor = FALSE)

# annotate data with protein, gene, and entrez ids
idx <- match(tidy_pep$Protein,gene_map$uniprot)
idy <- protein_identifier_col
tidy_pep <- tidy_pep %>% mutate(Entrez = gene_map$entrez[idx])
tidy_pep$Protein <- gene_map$uniprot[match(tidy_pep[[idy]],gene_map$id)]
tidy_pep$Entrez <- gene_map$entrez[match(tidy_pep[[idy]],gene_map$id)]
tidy_pep$Symbol <- gene_map$symbol[match(tidy_pep[[idy]],gene_map$id)]

# merge with sample meta data
tidy_pep <- tidy_pep %>% left_join(meta_data, by = "Sample")


## ---- insure that keratins have been removed

is_keratin <- grepl("keratin|Keratin", tidy_pep$"Protein Description")
keratins <- unique(tidy_pep$Protein[is_keratin])
tidy_pep <- tidy_pep %>% filter(!(Protein %in% keratins))

# Keratins removed: 22
message("Keratins removed: ", length(keratins))


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

# Mitochondrial proteins removed: 209
nMito <- sum(mito %in% tidy_pep$Entrez)
message("Mitochondrial proteins removed: ", nMito)

# rm mito prots
tidy_pep <- tidy_pep %>% dplyr::filter(!(Entrez %in% mito))
tidy_pep$Sample <- as.character(tidy_pep$Sample)
tidy_pep$Condition <- as.character(tidy_pep$Condition)


## ---- sample loading normalization 

# sl normalization
sl_pep <- normSL(tidy_pep)

# Following Sample loading normalization, 
# sample (run) Intensity sums are equal, e.g.:
# |Sample  | Total Intensity|
# |:-------|---------------:|
# |Colly-1 |    1.556415e+12|
# |Con-2   |    1.556415e+12|
# |Gephy-1 |    1.556415e+12|
# |SPQC-3  |    1.556415e+12|
# |UPF-1   |    1.556415e+12|

# Check, column sums (Run total intensity) should now be equal:
df <- sl_pep %>% group_by(Sample) %>% 
	summarize("Total Intensity"=sum(Intensity,na.rm=TRUE),.groups="drop")
knitr::kable(df)


## ---- protein summarization 

# protein summarization by sum
raw_prot <- sl_pep %>% group_by(Protein, Sample) %>% 
	summarize(Peptides = sum(!is.na(Intensity)), 
		  Intensity = sum(Intensity,na.rm=TRUE),.groups="drop") %>%
        left_join(meta_data,by="Sample")


## ---- protein level filtering 

# Remove one hit wonders -- proteins identified by a single peptide
one_hit_wonders <- unique(raw_prot$Protein[raw_prot$Peptides == 1])

# Number of one-hit-wonders: 674
message("Number of one-hit-wonders: ", length(one_hit_wonders))

# rm one hit wonders
filt_prot <- raw_prot %>% filter(!(Protein %in% one_hit_wonders))


## ---- sample pool normalization 
# perform normalization to QC samples

spn_prot <- normSP(filt_prot, pool="QC")


## ---- final normalized protein-level data

# insure 0 is NA 
spn_prot$Intensity[spn_prot$Intensity == 0] <- NA

# cast to dm and check for missing
dm <- spn_prot %>% dcast(Protein ~ Sample, value.var = "Intensity")
idx <- apply(dm,1,function(x) any(is.na(x)))
to_rm <- dm$Protein[idx]

message("rm proteins with incomplete cases: ", length(to_rm))

# remove the QC data before statistical testing:
tidy_prot <- spn_prot %>% filter(Condition != "QC") %>%
       	filter(!(Protein %in% to_rm))


## ---- statistical testing for positive control, e.g. Gephryin

# get uniprot ids
gphn <- geneLists::mapIDs("Gphn","symbol","uniprot",gene_map)
insyn1 <- geneLists::mapIDs("Insyn1","symbol","uniprot",gene_map)
collybistin <- geneLists::mapIDs("Arhgef9","symbol","uniprot",gene_map)

# check the final data:
knitr::kable(tidy_prot %>% subset(Protein == gphn))

# there are several ways we could go about setting up the models and contrasts:
# fx <- log2(Intensity) ~ 0 + Condition # intercept = 0
# fx <- log2(Intensity) ~ 1 + Condition # intercept = 1
# fx <- log2(Intensity) ~ Condition  # set Control as first level of Condition

# simple linear model:
fx <- log2(Intensity) ~ 0 + Condition

# NOTE: by setting the intercept to 0 we explicitly estimate all levels of
# Condition

# you may also fit:
#fx <- log2(Intensity) ~ Condition
#ordered_levels <- c("Control",
#		    "Collybistin-BioID",
#		    "Gephryin-BioID",
#		    "InSyn1-BioID")
#tidy_prot <- tidy_prot %>% 
#	mutate(Condition = factor(Condition, levels=ordered_levels))
#stopifnot(levels(tidy_prot$Condition)[1] == "Control")

# fit the model for gphn, insyn1, and collybistin
fm1 <- lm(fx, data = tidy_prot %>% subset(Protein == gphn))
fm2 <- lm(fx, data = tidy_prot %>% subset(Protein == insyn1))
fm3 <- lm(fx, data = tidy_prot %>% subset(Protein == collybistin))

# create Gephyrin vs BioID control contrast:
LT1 <- tidyProt::getContrast(fm1, "Gephyrin","Control")

# create Insyn1 contrast:
LT2 <- tidyProt::getContrast(fm2, "InSyn1","Control")

# create Arhgef9 contrast:
LT3 <- tidyProt::getContrast(fm3, "Collybistin","Control")

# assess contrast for gphn 
res1 <- lmTestContrast(fm1, LT1) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% 
	mutate(Symbol = geneLists::mapIDs(gphn,"uniprot","symbol",gene_map)) %>%
	dplyr::select(Symbol, everything())

# insyn1
res2 <- lmTestContrast(fm2, LT2) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% 
	mutate(Symbol=mapIDs(insyn1,"uniprot","symbol",gene_map)) %>%
	dplyr::select(Symbol, everything())

# arhgef9
res3 <- lmTestContrast(fm3, LT3) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% 
	mutate(Symbol = mapIDs(collybistin,"uniprot","symbol",gene_map)) %>%
	dplyr::select(Symbol, everything())

# check results
bind_rows(res1,res2,res3) %>% knitr::kable()



## ---- loop to perform test for all proteins

# empty lists for models, sigma2, and degrees of freedom for each fit
fm_list <- list()
s2_list <- list()
df_list <- list()

# fit the models: 
proteins <- unique(tidy_prot$Protein)
for (prot in proteins) {
	fm <- lm(fx, tidy_prot %>% subset(Protein == prot))
	fm_list[[prot]] <- fm
	s2_list[[prot]] <- sigma(fm)^2 # calculate sigma2 for moderation
	df_list[[prot]] <- fm$df.residual # get df for moderation
}

# perform t-statistic moderation
eb_var <- limma::squeezeVar(unlist(s2_list), unlist(df_list))
df_prior <- eb_var$df.prior
s2_prior <- eb_var$s2.prior
if (is.null(s2_prior)) { s2_prior <- 0 }

# s2 prior == 0, interpretation??

# examine results again with moderation
mod_res1 <- lmTestContrast(fm1, LT1, s2_prior, df_prior) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% 
	mutate(Symbol = mapIDs(gphn,"uniprot","symbol",gene_map)) %>%
	dplyr::select(Symbol, everything())
mod_res2 <- lmTestContrast(fm2, LT2, s2_prior, df_prior) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% 
	mutate(Symbol = mapIDs(insyn1,"uniprot","symbol",gene_map)) %>%
	dplyr::select(Symbol, everything())
mod_res3 <- lmTestContrast(fm3, LT3, s2_prior, df_prior) %>% 
	mutate(Pvalue=formatC(Pvalue)) %>% 
	mutate(Symbol = mapIDs(collybistin,"uniprot","symbol",gene_map)) %>%
	dplyr::select(Symbol, everything())

# check moderated results
bind_rows(mod_res1,mod_res2,mod_res3) %>% 
	knitr::kable()

# loop to moderate comparisons for all proteins
result_list <- list()
for (prot in proteins) {
	fm <- fm_list[[prot]]
	# assess all three types of contrast, pass s2_prior and df_prior 
	# to moderate the comparison
	result_list[[prot]] <- do.call(rbind, 
				       lapply(list(LT1,LT2,LT3), function(x) 
		                       lmTestContrast(fm,x,s2_prior,df_prior)))
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
idx <- match(results_df$Protein,gene_map$uniprot)
results_df <- results_df %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx], .after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx], .after="Symbol")

# final sort
results_df <- results_df %>% 
	arrange(desc(sig), desc(up), Pvalue, desc(log2FC)) %>%
	dplyr::select(-sig,-up)

# total instances of sig enrichment:
data.table("nSig"=sum(results_df$candidate)) %>% knitr::kable()

# summary for each contrast
results_df %>% 
	group_by(Contrast) %>% 
	summarize(nSig=sum(candidate),.groups="drop") %>% knitr::kable()

# overlap
sig_prots <- split(results_df$Protein[results_df$candidate],
		   results_df$Contrast[results_df$candidate])
overlap <- Reduce(intersect,sig_prots)
data.table("overlap" = length(overlap)) %>% knitr::kable()

# common prots:
results_df %>% 
	filter(Protein %in% overlap) %>% 
	arrange(Symbol,Contrast) %>% 
	knitr::kable()


# add to list
sig_prots[["overlap"]] <- overlap


## ---- save results

# annot
idx <- match(spn_prot$Protein,gene_map$uniprot)
spn_prot$Symbol <- gene_map$symbol[idx]
spn_prot$Entrez <- gene_map$entrez[idx]

# split into list
results_list <- results_df %>% group_by(Contrast) %>% group_split()

# simplify list names --> sheet names
namen <- gsub("Condition|Control|-","",
	      sapply(results_list, function(x) unique(x$Contrast)))
names(results_list) <- namen

class(results_list) <- "list"

# final results_list to be saved as excel document
results_list[["Raw Data"]] <- as.data.frame(raw_pep)
results_list[["Norm Data"]] <- spn_prot

# save tidy, normalized protein data as rda
myfile <- file.path(root,"data","ipsd_bioid.rda") 
ipsd_bioid <- spn_prot
save(ipsd_bioid,file=myfile,version=2)
message("saved: ", myfile)

# save statistical results as excel
myfile <- file.path(root,"tables","iPSD_iBioID_Results.xlsx")
write_excel(results_list, myfile)
message("saved: ", myfile)

# save gene_map
myfile <- file.path(root,"data","ipsd_gene_map.rda")
save(gene_map,file=myfile,version=2)
message("saved: ", myfile)

# save statistical results 
myfile <- file.path(root,"data","ipsd_results.rda")
ipsd_results <- results_df
save(ipsd_results,file=myfile,version=2)
message("saved: ", myfile)

# save sig_prots
myfile <- file.path(root,"data","ipsd_sig_prots.rda")
save(sig_prots,file=myfile,version=2)
message("saved: ", myfile)

# save key genes 
myfile <- file.path(root,"data","gphn.rda")
save(gphn,file=myfile,version=2)
message("saved: ", myfile)
#
myfile <- file.path(root,"data","insyn1.rda")
save(insyn1,file=myfile,version=2)
message("saved: ", myfile)
#
myfile <- file.path(root,"data","collybistin.rda")
save(collybistin,file=myfile,version=2)
message("saved: ", myfile)
#
myfile <- file.path(root,"data","insyn2.rda")
save(insyn2,file=myfile,version=2)
message("saved: ", myfile)
