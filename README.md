# Uezu2016

iBioID proteoimcs data from Uezu _et al.,_ 2016, published in
[Science](https://science.sciencemag.org/content/353/6304/1123.full).

To access the raw data in R, install the package with devtools, and then access
the data in `inst/extdta` using `system.file` and
[readxl](https://readxl.tidyverse.org/) (see example below). Or, if that sounds
complicated, you can simply download the raw files from this repository,
[here](./inst/extdata).

```R

devtools::install_github("soderling-lab/Uezu2016")

library(readxl)

# the Uezu2016 package contains the following files in inst/extdata:
epsd <- "ePSD_BioID-Supplementary-data.xlsx"
ipsd <- "iPSD_BioID-Supplementary-data.xlsx"

# access your system's path to the data
myfile <- system.file("extdata", ipsd, package="Uezu2016")

readxl::excel_sheets(myfile)

# [1] "Table 1" "Table 2" "Table 3" 
# ^MS info, Sample info, and the raw (peptide-level) data

df <- readxl::read_excel(myfile, sheet = 3) 

```

![ipsd](./ipsd.png)


## Citation

```bibtex
@article{Uezu2016,
	doi = {10.1126/science.aag0821},
	url = {https://doi.org/10.1126%2Fscience.aag0821},
	year = 2016,
	month = {sep},
	publisher = {American Association for the Advancement of Science ({AAAS})},
	volume = {353},
	number = {6304},
	pages = {1123--1129},
	author = {A. Uezu and D. J. Kanak and T. W. A. Bradshaw and E. J. Soderblom and C. M. Catavero and A. C. Burette and R. J. Weinberg and S. H. Soderling},
	title = {Identification of an elaborate complex mediating postsynaptic inhibition},
	journal = {Science}
}
```
