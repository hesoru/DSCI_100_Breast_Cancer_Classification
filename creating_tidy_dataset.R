setwd("/Users/Helena/Documents/BCS 1ST YEAR/DSCI 100/Group Project/breast_cancer_proteome")

install.packages("/Users/Helena/Downloads/unpivotr_0.6.3.tar.gz", repos = NULL, type = "source")
install.packages("/Users/Helena/Downloads/janitor_2.2.0.tar.gz", repos = NULL, type = "source")

library("tidyverse")
library("stringr")
library("janitor")


# reading the proteome and clinical data CSVs, cleaning column names
proteome_data <- read_csv("77_cancer_proteomes_CPTAC_itraq.csv")

clinical_data <- read_csv("clinical_data_breast_cancer.csv") 
clinical_data <- clean_names(clinical_data)


# clinical data: change complete_tcga_id to uniform format
clinical_data[c("tcga", "tcga_code_1", "tcga_code_2")] <- str_split_fixed(clinical_data$complete_tcga_id, "-", 3)
clinical_data$tcga_id <- paste(clinical_data$tcga_code_1, clinical_data$tcga_code_2, sep = "-")
clinical_data_tidy <- clinical_data[-c(1, 31:33)] %>% relocate(tcga_id, .before = gender)


# proteome data: convert the donor columns into a single donor column
proteome_data_longer <- pivot_longer(proteome_data, 4:86, names_to = "complete_tcga_id", values_to = "protein_expression_log2_iTRAQ_ratios")
# proteome data: change complete_tcga_id to uniform format
proteome_data_longer[c("tcga_id", "tcga_id_2")] <- str_split_fixed(proteome_data_longer$complete_tcga_id, "\\.", 2)
proteome_data_longer_tidy <- proteome_data_longer[-c(4, 7)] %>% relocate(tcga_id, .before = RefSeq_accession_number)


# merge proteome and clinical data
proteome_and_clinical_data_tidy <- merge(proteome_data_longer_tidy, clinical_data_tidy, by="tcga_id")

# validate that there are 77 donors and 12,553 genes per donor
donors_and_genes_per_donor <- proteome_and_clinical_data_tidy %>% group_by(tcga_id) %>% count()
# 77 donors (tcga_id) and 12,553 genes (n) per donor

write.csv(proteome_and_clinical_data_tidy, "proteome_and_clinical_data_tidy.csv")
