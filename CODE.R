rm(list = ls())  
library(dplyr)
library(MRPRESSO)
library(TwoSampleMR)
library(dplyr)
library(ggsci)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggview)
library(plotly)
library(TSMRhelper)
library(forestploter)
library(grid)
library(stringr)
library(grDevices)
library(parallel)
library(vroom)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

# Read exposure data
exp_clumped <- vroom("338_all_data.csv")

# Extract outcome data
outcome_dat <- extract_outcome_data(snps = exp_clumped$SNP, outcomes = "ieu-b-7", proxies = FALSE)

# Harmonize data
dat <- harmonise_data(exposure_dat = exp_clumped, outcome_dat = outcome_dat)
dat <- distinct(dat, SNP, .keep_all = TRUE)
dat <- dat[dat$mr_keep,]
write.csv(dat, file = "338_PD_dat.csv", row.names = FALSE)

# Split data by exposure id
dat <- split(dat, dat$id.exposure)

# Run MR analysis
cl <- makeCluster(detectCores())
res <- parLapply(cl, dat, mr)
stopCluster(cl)
res <- do.call(rbind, res)
res <- generate_odds_ratios(res)
write.csv(res, file = "338_PD_res.csv", row.names = FALSE)

# Heterogeneity test
dat <- do.call(rbind, dat)
heterogeneity <- mr_heterogeneity(dat)
write.csv(heterogeneity, file = "338_PD_heterogeneity.csv", row.names = FALSE)

# Pleiotropy test
pleiotropy_test <- mr_pleiotropy_test(dat)
write.csv(pleiotropy_test, file = "338_PD_pleiotropy_test.csv", row.names = FALSE)

####### Reverse MR ######
exp_clumped <- extract_instruments("ieu-b-7", p1 = 5e-8)
table(!((((exp_clumped$beta.exposure)^2) / ((exp_clumped$se.exposure)^2)) < 10))
exp_clumped <- exp_clumped[!((((exp_clumped$beta.exposure)^2) / ((exp_clumped$se.exposure)^2)) < 10),]

# Read trait file
idTotrait <- read.csv("./338_Trait.csv")

# Read result id file
total_data <- read.csv("id_338_PD.csv")
id_list <- total_data$id

# Set folder path for raw data
folder_path <- "F:/SNPdata"

# Construct file paths and add suffix
file_paths <- file.path(folder_path, paste0(id_list, "_buildGRCh37.tsv.csv"))

dat <- list()

for (i in 1:length(id_list)) {
  print(id_list[i])
  file_path <- file_paths[i]
  outcome_dat <- vroom(file_path)
  outcome_dat$p_value <- as.numeric(outcome_dat$p_value)
  
  # Filter SNPs
  outcome_dat <- outcome_dat %>% filter(SNP %in% unique(exp_clumped$SNP))
  outcome_dat$id <- str_split(id_list[i], "_buildGRCh37.tsv.csv", simplify = TRUE)[, 1]
  outcome_dat$phe <- idTotrait[idTotrait$id == outcome_dat$id[1], ]$Trait
  
  # Format data
  outcome_dat <- format_data(outcome_dat,
                             type = "outcome",
                             snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "standard_error",
                             eaf_col = "effect_allele_frequency",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "p_value",
                             chr_col = "CHR",
                             pos_col = "BP",
                             id_col = "id",
                             phenotype_col = "phe")
  dat[[i]] <- harmonise_data(exposure_dat = exp_clumped, outcome_dat = outcome_dat)
}
dat <- do.call(rbind, dat)

dat <- dat[dat$mr_keep,]
write.csv(dat, file = "PD_338_dat.csv", row.names = FALSE)

re_res <- mr(dat)
re_res <- generate_odds_ratios(re_res)
write.csv(re_res, file = "PD_338_res.csv", row.names = FALSE)

# Heterogeneity test
heterogeneity <- mr_heterogeneity(dat)
write.csv(heterogeneity, file = "PD_338_heterogeneity.csv", row.names = FALSE)

# Pleiotropy test
pleiotropy_test <- mr_pleiotropy_test(dat)
write.csv(pleiotropy_test, file = "PD_338_pleiotropy_test.csv", row.names = FALSE)
