# Uezu2016

iBioID proteoimcs data from Uezu __et al.,__ 2016, published in
[Science](https://science.sciencemag.org/content/353/6304/1123.full).

Access the raw data in R using [readxl](https://readxl.tidyverse.org/).
```R
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
