# Tutorial for GWAS Data Liftover (from hg19 to hg38)
*Created by:* Yazdan Asgari<br>
*Creation date:* 24 Nov 2023<br>
*Update:* Feb 2024<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>

This repository contains an R script that performs the liftover of GWAS (Genome-Wide Association Study) data from hg19 to hg38 coordinates using the `rtracklayer`, `GenomicRanges`, and `vroom` packages. The liftover is achieved through the use of a chain file provided by UCSC.

## Prerequisites
Before running the script, make sure you have the required R packages installed:

```R
install.packages(c("rtracklayer", "GenomicRanges", "vroom"))
```

**IMPORTANT NOTE:** This script needs these two column names available in your GWAS data (*"chr"* and *"bp_hg19"*). *"chr"* column is in **number format**. If their names are different in your input data, do not forget to change them in the script before using it.

## Data Processing Steps
1. Reading GWAS data
```R
library(rtracklayer)
library(GenomicRanges)
library(vroom)

path_gwas <- "/path/to/your/GWAS/"
gwas_hg19 <- vroom(file = paste0(path_gwas, "GWAS_hg19.txt"))
head(gwas_hg19)
gwas_hg19['chr_str'] <- paste0("chr", gwas_hg19$chr)

gwas_hg19_sel <- gwas_hg19

# IMPORTANT NOTE: We saw that sometimes if you did not add "strand" column info on your data, some bp_hg38 positions will add with +1 bp in the created hg38 file !!! So, we added this column to the data as well.
gwas_hg19_sel["strand"] <- "+"

# Check data types of gwas_hg19_sel
str(gwas_hg19_sel)

# Use this section ONLY if you have "NA" data in "bp_hg19" column
# Check for missing values in "bp_hg19" column. If so, it then removes them
#sum(is.na(gwas_hg19_sel$bp_hg19))
#gwas_hg19_sel <- gwas_hg19_sel[!is.na(gwas_hg19_sel$bp_hg19),]
#sum(is.na(gwas_hg19_sel$bp_hg19))

```

2. Importing chain file
```R
# It was downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/
chain <- import.chain("hg19ToHg38.over.chain")
```
You could also download the "chain file" from the [0_files](/0_files) folder.
<br>
<br>
3. Creating GRanges object
```R
gr <- makeGRangesFromDataFrame(gwas_hg19_sel, ignore.strand = TRUE, seqnames.field = "chr_str", start.field = "bp_hg19", end.field = "bp_hg19")
```

4. Lifting over to hg38
```R
hg38 <- liftOver(gr, chain)
hg38_df <- as.data.frame(hg38)
```

5. Matching and updating data
```R
gwas_hg19_sel$rownames_col <- as.numeric(rownames(gwas_hg19_sel))
head(gwas_hg19_sel)
match_pos <- match(gwas_hg19_sel$rownames_col, hg38_df$group)

# Add the corresponding information from hg38_df to gwas_hg19_sel
gwas_hg19_sel$group[!is.na(match_pos)] <- hg38_df$group[match_pos[!is.na(match_pos)]]
gwas_hg19_sel$seqnames[!is.na(match_pos)] <- hg38_df$seqnames[match_pos[!is.na(match_pos)]]
gwas_hg19_sel$bp_hg38[!is.na(match_pos)] <- hg38_df$start[match_pos[!is.na(match_pos)]]
head(gwas_hg19_sel)
```

6. Checking for missing values and writing output
```R
# Check for missing values in the "bp_hg38" column of "gwas_hg19_sel" GWAS file, then remove them
sum(is.na(gwas_hg19_sel$bp_hg38))
gwas_hg19_sel <- gwas_hg19_sel[!is.na(gwas_hg19_sel$bp_hg38),]
sum(is.na(gwas_hg19_sel$bp_hg38))
# Write the output to a file
vroom_write(gwas_hg19_sel, file = "GWAS_hg38.txt", delim = "\t", quote = "none", col_names = TRUE)
```

**Notes**
- Ensure that the hg19 to hg38 chain file is present in the project directory.
- Verify that the required R packages are installed before running the script.

# Tutorial for GWAS Data Liftover (from hg38 to hg19)
Just simply use the exact previous procedure, but this time use *"hg38ToHg19.over.chain"* chain file (available in the [0_files](/0_files) folder) and switch the variable names *"hg19"* and *"hg38"* in the script.
