library(tidyverse)

########################################################
# Generating Ecoli dataframe 
########################################################

################ Dorado 5mC dataframes #################

Ecoli_MSssI_5mC_bisulfite_r10 <- read.csv(file = "bisulfite/meta_20x/Ecoli_DM_MSssI_rep01_20x.bed",sep = '\t', header = T)
Ecoli_MCviP_5mC_bisulfite_r10 <- read.csv(file = "bisulfite/meta_20x/Ecoli_DM_MCviP_rep01_20x.bed",sep = '\t', header = T)
Ecoli_DM_5mC_bisulfite_r10 <- read.csv(file = "bisulfite/meta_20x/Ecoli_DM_rep01_20x.bed",sep = '\t', header = T)
Ecoli_WT_5mC_bisulfite_r10 <- read.csv(file = "bisulfite/meta_20x/Ecoli_WT_rep01_20x.bed",sep = '\t', header = T)

Ecoli_MSssI_5mC_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_4kHz_sup_4kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_MCviP_5mC_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MCviP_rep01_kit14_4kHz_sup_4kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_DM_5mC_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_4kHz_sup_4kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_WT_5mC_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_4kHz_sup_4kHz_v4_5mC.bed",sep = '\t', header = T)

Ecoli_MSssI_5mC_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_sup_5kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_MCviP_5mC_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MCviP_rep01_kit14_5kHz_sup_5kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_DM_5mC_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_sup_5kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_WT_5mC_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_sup_5kHz_v4_5mC.bed",sep = '\t', header = T)

Ecoli_MSssI_5mC_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_hac_5kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_MCviP_5mC_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MCviP_rep01_kit14_5kHz_hac_5kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_DM_5mC_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_hac_5kHz_v4_5mC.bed",sep = '\t', header = T)
Ecoli_WT_5mC_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_hac_5kHz_v4_5mC.bed",sep = '\t', header = T)


Ecoli_MSssI_5mC_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_sup_5kHz_v5_5mC.bed",sep = '\t', header = T)
Ecoli_MCviP_5mC_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MCviP_rep01_kit14_5kHz_sup_5kHz_v5_5mC.bed",sep = '\t', header = T)
Ecoli_DM_5mC_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_sup_5kHz_v5_5mC.bed",sep = '\t', header = T)
Ecoli_WT_5mC_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_sup_5kHz_v5_5mC.bed",sep = '\t', header = T)

Ecoli_MSssI_5mC_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_hac_5kHz_v5_5mC.bed",sep = '\t', header = T)
Ecoli_MCviP_5mC_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MCviP_rep01_kit14_5kHz_hac_5kHz_v5_5mC.bed",sep = '\t', header = T)
Ecoli_DM_5mC_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_hac_5kHz_v5_5mC.bed",sep = '\t', header = T)
Ecoli_WT_5mC_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_hac_5kHz_v5_5mC.bed",sep = '\t', header = T)


Ecoli_WT_5mC_bisulfite_r10 <- Ecoli_WT_5mC_bisulfite_r10 %>% mutate(sample="Ecoli_WT")
Ecoli_MSssI_5mC_dorado_4khz_r10 <- Ecoli_MSssI_5mC_dorado_4khz_r10 %>% mutate(sample="Ecoli_MSssI")
Ecoli_MCviP_5mC_dorado_4khz_r10 <- Ecoli_MCviP_5mC_dorado_4khz_r10 %>% mutate(sample="Ecoli_MCviP")
Ecoli_DM_5mC_dorado_4khz_r10 <- Ecoli_DM_5mC_dorado_4khz_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_5mC_dorado_4khz_r10 <- Ecoli_WT_5mC_dorado_4khz_r10 %>% mutate(sample="Ecoli_WT")

Ecoli_MSssI_5mC_dorado_5khz_v4sup_r10 <- Ecoli_MSssI_5mC_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_MSssI")
Ecoli_MCviP_5mC_dorado_5khz_v4sup_r10 <- Ecoli_MCviP_5mC_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_MCviP")
Ecoli_DM_5mC_dorado_5khz_v4sup_r10 <- Ecoli_DM_5mC_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_5mC_dorado_5khz_v4sup_r10 <- Ecoli_WT_5mC_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_WT")

Ecoli_MSssI_5mC_dorado_5khz_v4hac_r10 <- Ecoli_MSssI_5mC_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_MSssI")
Ecoli_MCviP_5mC_dorado_5khz_v4hac_r10 <- Ecoli_MCviP_5mC_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_MCviP")
Ecoli_DM_5mC_dorado_5khz_v4hac_r10 <- Ecoli_DM_5mC_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_5mC_dorado_5khz_v4hac_r10 <- Ecoli_WT_5mC_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_WT")

Ecoli_MSssI_5mC_dorado_5khz_v5sup_r10 <- Ecoli_MSssI_5mC_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_MSssI")
Ecoli_MCviP_5mC_dorado_5khz_v5sup_r10 <- Ecoli_MCviP_5mC_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_MCviP")
Ecoli_DM_5mC_dorado_5khz_v5sup_r10 <- Ecoli_DM_5mC_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_5mC_dorado_5khz_v5sup_r10 <- Ecoli_WT_5mC_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_WT")

Ecoli_MSssI_5mC_dorado_5khz_v5hac_r10 <- Ecoli_MSssI_5mC_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_MSssI")
Ecoli_MCviP_5mC_dorado_5khz_v5hac_r10 <- Ecoli_MCviP_5mC_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_MCviP")
Ecoli_DM_5mC_dorado_5khz_v5hac_r10 <- Ecoli_DM_5mC_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_5mC_dorado_5khz_v5hac_r10 <- Ecoli_WT_5mC_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_WT")


Ecoli_5mC_r10_meta <- rbind(Ecoli_MSssI_5mC_bisulfite_r10,Ecoli_MCviP_5mC_bisulfite_r10,Ecoli_DM_5mC_bisulfite_r10,Ecoli_WT_5mC_bisulfite_r10,
                            Ecoli_MSssI_5mC_dorado_4khz_r10,Ecoli_MCviP_5mC_dorado_4khz_r10,Ecoli_DM_5mC_dorado_4khz_r10,Ecoli_WT_5mC_dorado_4khz_r10,
                            Ecoli_MSssI_5mC_dorado_5khz_v4sup_r10,Ecoli_MCviP_5mC_dorado_5khz_v4sup_r10,Ecoli_DM_5mC_dorado_5khz_v4sup_r10,Ecoli_WT_5mC_dorado_5khz_v4sup_r10,
                            Ecoli_MSssI_5mC_dorado_5khz_v4hac_r10,Ecoli_MCviP_5mC_dorado_5khz_v4hac_r10,Ecoli_DM_5mC_dorado_5khz_v4hac_r10,Ecoli_WT_5mC_dorado_5khz_v4hac_r10,
                            Ecoli_MSssI_5mC_dorado_5khz_v5sup_r10,Ecoli_MCviP_5mC_dorado_5khz_v5sup_r10,Ecoli_DM_5mC_dorado_5khz_v5sup_r10,Ecoli_WT_5mC_dorado_5khz_v5sup_r10,
                            Ecoli_MSssI_5mC_dorado_5khz_v5hac_r10,Ecoli_MCviP_5mC_dorado_5khz_v5hac_r10,Ecoli_DM_5mC_dorado_5khz_v5hac_r10,Ecoli_WT_5mC_dorado_5khz_v5hac_r10) 


################ Dorado 6mA dataframes #################

Ecoli_DM_6mA_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_4kHz_sup_4kHz_v4_6mA.bed",sep = '\t', header = T)
Ecoli_WT_6mA_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_4kHz_sup_4kHz_v4_6mA.bed",sep = '\t', header = T)
Ecoli_MSssI_6mA_dorado_4khz_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_4kHz_sup_4kHz_v4_6mA.bed",sep = '\t', header = T)


Ecoli_DM_6mA_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_sup_5kHz_v4_6mA.bed",sep = '\t', header = T)
Ecoli_WT_6mA_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_sup_5kHz_v4_6mA.bed",sep = '\t', header = T)
Ecoli_MSssI_6mA_dorado_5khz_v4sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_sup_5kHz_v4_6mA.bed",sep = '\t', header = T)


Ecoli_DM_6mA_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_hac_5kHz_v4_6mA.bed",sep = '\t', header = T)
Ecoli_WT_6mA_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_hac_5kHz_v4_6mA.bed",sep = '\t', header = T)
Ecoli_MSssI_6mA_dorado_5khz_v4hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_hac_5kHz_v4_6mA.bed",sep = '\t', header = T)


Ecoli_DM_6mA_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_sup_5kHz_v5_6mA.bed",sep = '\t', header = T)
Ecoli_WT_6mA_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_sup_5kHz_v5_6mA.bed",sep = '\t', header = T)
Ecoli_MSssI_6mA_dorado_5khz_v5sup_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_sup_5kHz_v5_6mA.bed",sep = '\t', header = T)


Ecoli_DM_6mA_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_no_rep01_kit14_5kHz_hac_5kHz_v5_6mA.bed",sep = '\t', header = T)
Ecoli_WT_6mA_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_WT_NAT_no_rep01_kit14_5kHz_hac_5kHz_v5_6mA.bed",sep = '\t', header = T)
Ecoli_MSssI_6mA_dorado_5khz_v5hac_r10 <- read.csv(file = "Ecoli_meta_20x/Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_hac_5kHz_v5_6mA.bed",sep = '\t', header = T)


Ecoli_DM_6mA_dorado_4khz_r10 <- Ecoli_DM_6mA_dorado_4khz_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_6mA_dorado_4khz_r10 <- Ecoli_WT_6mA_dorado_4khz_r10 %>% mutate(sample="Ecoli_WT")
Ecoli_MSssI_6mA_dorado_4khz_r10 <- Ecoli_MSssI_6mA_dorado_4khz_r10 %>% mutate(sample="Ecoli_MSssI")

Ecoli_DM_6mA_dorado_5khz_v4sup_r10 <- Ecoli_DM_6mA_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_6mA_dorado_5khz_v4sup_r10 <- Ecoli_WT_6mA_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_WT")
Ecoli_MSssI_6mA_dorado_5khz_v4sup_r10 <- Ecoli_MSssI_6mA_dorado_5khz_v4sup_r10 %>% mutate(sample="Ecoli_MSssI")

Ecoli_DM_6mA_dorado_5khz_v4hac_r10 <- Ecoli_DM_6mA_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_6mA_dorado_5khz_v4hac_r10 <- Ecoli_WT_6mA_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_WT")
Ecoli_MSssI_6mA_dorado_5khz_v4hac_r10 <- Ecoli_MSssI_6mA_dorado_5khz_v4hac_r10 %>% mutate(sample="Ecoli_MSssI")

Ecoli_DM_6mA_dorado_5khz_v5sup_r10 <- Ecoli_DM_6mA_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_6mA_dorado_5khz_v5sup_r10 <- Ecoli_WT_6mA_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_WT")
Ecoli_MSssI_6mA_dorado_5khz_v5sup_r10 <- Ecoli_MSssI_6mA_dorado_5khz_v5sup_r10 %>% mutate(sample="Ecoli_MSssI")

Ecoli_DM_6mA_dorado_5khz_v5hac_r10 <- Ecoli_DM_6mA_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_DM")
Ecoli_WT_6mA_dorado_5khz_v5hac_r10 <- Ecoli_WT_6mA_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_WT")
Ecoli_MSssI_6mA_dorado_5khz_v5hac_r10 <- Ecoli_MSssI_6mA_dorado_5khz_v5hac_r10 %>% mutate(sample="Ecoli_MSssI")


Ecoli_6mA_r10_meta <- rbind(Ecoli_DM_6mA_dorado_4khz_r10,Ecoli_WT_6mA_dorado_4khz_r10,Ecoli_MSssI_6mA_dorado_4khz_r10,
                            Ecoli_DM_6mA_dorado_5khz_v4sup_r10,Ecoli_WT_6mA_dorado_5khz_v4sup_r10,Ecoli_MSssI_6mA_dorado_5khz_v4sup_r10,
                            Ecoli_DM_6mA_dorado_5khz_v4hac_r10,Ecoli_WT_6mA_dorado_5khz_v4hac_r10,Ecoli_MSssI_6mA_dorado_5khz_v4hac_r10,
                            Ecoli_DM_6mA_dorado_5khz_v5sup_r10,Ecoli_WT_6mA_dorado_5khz_v5sup_r10,Ecoli_MSssI_6mA_dorado_5khz_v5sup_r10,
                            Ecoli_DM_6mA_dorado_5khz_v5hac_r10,Ecoli_WT_6mA_dorado_5khz_v5hac_r10,Ecoli_MSssI_6mA_dorado_5khz_v5hac_r10) 


################ Deepmod 5mC dataframes #################


Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_bilstm <- read.csv(file = "Ecoli_DM_NAT_MSssI_rep01_kit14_4kHz_100x_bilstmv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_bilstm <- read.csv(file = "Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_100x_bilstmv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_DM_NAT_no_deepmod_4khz_r10_bilstm <- read.csv(file = "Ecoli_DM_NAT_no_rep01_kit14_4kHz_100x_bilstmv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_DM_NAT_no_deepmod_5khz_r10_bilstm <- read.csv(file = "Ecoli_DM_NAT_no_rep01_kit14_5kHz_100x_bilstmv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_WT_NAT_no_deepmod_4khz_r10_bilstm <- read.csv(file = "Ecoli_WT_NAT_no_rep01_kit14_4kHz_100x_bilstmv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_WT_NAT_no_deepmod_5khz_r10_bilstm <- read.csv(file = "Ecoli_WT_NAT_no_rep01_kit14_5kHz_100x_bilstmv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)


Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_transformer <- read.csv(file = "Ecoli_DM_NAT_MSssI_rep01_kit14_4kHz_100x_transformerv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_transformer <- read.csv(file = "Ecoli_DM_NAT_MSssI_rep01_kit14_5kHz_100x_transformerv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_DM_NAT_no_deepmod_4khz_r10_transformer <- read.csv(file = "Ecoli_DM_NAT_no_rep01_kit14_4kHz_100x_transformerv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_DM_NAT_no_deepmod_5khz_r10_transformer <- read.csv(file = "Ecoli_DM_NAT_no_rep01_kit14_5kHz_100x_transformerv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_WT_NAT_no_deepmod_4khz_r10_transformer <- read.csv(file = "Ecoli_WT_NAT_no_rep01_kit14_4kHz_100x_transformerv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
Ecoli_WT_NAT_no_deepmod_5khz_r10_transformer <- read.csv(file = "Ecoli_WT_NAT_no_rep01_kit14_5kHz_100x_transformerv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)


Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_bilstm <- Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="4kHz",sample="Ecoli_MSssI",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_bilstm <- Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="5kHz",sample="Ecoli_MSssI",device="grid",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_DM_NAT_no_deepmod_4khz_r10_bilstm <- Ecoli_DM_NAT_no_deepmod_4khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="4kHz",sample="Ecoli_DM",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_DM_NAT_no_deepmod_5khz_r10_bilstm <- Ecoli_DM_NAT_no_deepmod_5khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="5kHz",sample="Ecoli_DM",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_WT_NAT_no_deepmod_4khz_r10_bilstm <- Ecoli_WT_NAT_no_deepmod_4khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="4kHz",sample="Ecoli_WT",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_WT_NAT_no_deepmod_5khz_r10_bilstm <- Ecoli_WT_NAT_no_deepmod_5khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="5kHz",sample="Ecoli_WT",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")


Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_transformer <- Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="4kHz",sample="Ecoli_MSssI",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_transformer <- Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="5kHz",sample="Ecoli_MSssI",device="grid",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_DM_NAT_no_deepmod_4khz_r10_transformer <- Ecoli_DM_NAT_no_deepmod_4khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="4kHz",sample="Ecoli_DM",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_DM_NAT_no_deepmod_5khz_r10_transformer <- Ecoli_DM_NAT_no_deepmod_5khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="5kHz",sample="Ecoli_DM",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_WT_NAT_no_deepmod_4khz_r10_transformer <- Ecoli_WT_NAT_no_deepmod_4khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="4kHz",sample="Ecoli_WT",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

Ecoli_WT_NAT_no_deepmod_5khz_r10_transformer <- Ecoli_WT_NAT_no_deepmod_5khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="5kHz",sample="Ecoli_WT",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="Ecoli") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")


Ecoli_deepmod_r10 <- rbind(Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_bilstm,Ecoli_DM_NAT_no_deepmod_4khz_r10_bilstm,Ecoli_WT_NAT_no_deepmod_4khz_r10_bilstm,
                          Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_bilstm,Ecoli_DM_NAT_no_deepmod_5khz_r10_bilstm,Ecoli_WT_NAT_no_deepmod_5khz_r10_bilstm,
                          Ecoli_DM_NAT_MSssI_deepmod_4khz_r10_transformer,Ecoli_DM_NAT_no_deepmod_4khz_r10_transformer,Ecoli_WT_NAT_no_deepmod_4khz_r10_transformer,
                          Ecoli_DM_NAT_MSssI_deepmod_5khz_r10_transformer,Ecoli_DM_NAT_no_deepmod_5khz_r10_transformer,Ecoli_WT_NAT_no_deepmod_5khz_r10_transformer)


Ecoli_allTools <- rbind(Ecoli_5mC_r10_meta,Ecoli_6mA_r10_meta,Ecoli_deepmod_r10)

write.table(Ecoli_allTools, file = "Ecoli_allTools_meta.tsv",sep ="\t",row.names=FALSE,quote = F)


########################################################
# Generating HP26695 dataframe 
########################################################

################ Dorado 5mC dataframes #################


HP26695_NAT_5mC_dorado_4khz_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_4kHz_sup_4kHz_v4_5mC.bed",sep = '\t', header = T)
HP26695_WGA_5mC_dorado_4khz_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_4kHz_sup_4kHz_v4_5mC.bed",sep = '\t', header = T)

HP26695_NAT_5mC_dorado_5khz_v4sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_sup_5kHz_v4_5mC.bed",sep = '\t', header = T)
HP26695_WGA_5mC_dorado_5khz_v4sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_sup_5kHz_v4_5mC.bed",sep = '\t', header = T)

HP26695_NAT_5mC_dorado_5khz_v4hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_hac_5kHz_v4_5mC.bed",sep = '\t', header = T)
HP26695_WGA_5mC_dorado_5khz_v4hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_hac_5kHz_v4_5mC.bed",sep = '\t', header = T)

HP26695_NAT_5mC_dorado_5khz_v5sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_sup_5kHz_v5_5mC.bed",sep = '\t', header = T)
HP26695_WGA_5mC_dorado_5khz_v5sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_sup_5kHz_v5_5mC.bed",sep = '\t', header = T)

HP26695_NAT_5mC_dorado_5khz_v5hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_hac_5kHz_v5_5mC.bed",sep = '\t', header = T)
HP26695_WGA_5mC_dorado_5khz_v5hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_hac_5kHz_v5_5mC.bed",sep = '\t', header = T)


HP26695_NAT_5mC_dorado_4khz_r10 <- HP26695_NAT_5mC_dorado_4khz_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_5mC_dorado_4khz_r10 <- HP26695_WGA_5mC_dorado_4khz_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_5mC_dorado_5khz_v4sup_r10 <- HP26695_NAT_5mC_dorado_5khz_v4sup_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_5mC_dorado_5khz_v4sup_r10 <- HP26695_WGA_5mC_dorado_5khz_v4sup_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_5mC_dorado_5khz_v4hac_r10 <- HP26695_NAT_5mC_dorado_5khz_v4hac_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_5mC_dorado_5khz_v4hac_r10 <- HP26695_WGA_5mC_dorado_5khz_v4hac_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_5mC_dorado_5khz_v5sup_r10 <- HP26695_NAT_5mC_dorado_5khz_v5sup_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_5mC_dorado_5khz_v5sup_r10 <- HP26695_WGA_5mC_dorado_5khz_v5sup_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_5mC_dorado_5khz_v5hac_r10 <- HP26695_NAT_5mC_dorado_5khz_v5hac_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_5mC_dorado_5khz_v5hac_r10 <- HP26695_WGA_5mC_dorado_5khz_v5hac_r10 %>% mutate(sample="HP26695_WGA")


HP26695_5mC_r10_meta <- rbind(HP26695_NAT_5mC_dorado_4khz_r10,HP26695_WGA_5mC_dorado_4khz_r10,
                              HP26695_NAT_5mC_dorado_5khz_v4sup_r10,HP26695_WGA_5mC_dorado_5khz_v4sup_r10,
                              HP26695_NAT_5mC_dorado_5khz_v4hac_r10,HP26695_WGA_5mC_dorado_5khz_v4hac_r10,
                              HP26695_NAT_5mC_dorado_5khz_v5sup_r10,HP26695_WGA_5mC_dorado_5khz_v5sup_r10,
                              HP26695_NAT_5mC_dorado_5khz_v5hac_r10,HP26695_WGA_5mC_dorado_5khz_v5hac_r10)


################ Dorado 6mA dataframes #################

HP26695_NAT_6mA_dorado_4khz_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_4kHz_sup_4kHz_v4_6mA.bed",sep = '\t', header = T)
HP26695_WGA_6mA_dorado_4khz_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_4kHz_sup_4kHz_v4_6mA.bed",sep = '\t', header = T)

HP26695_NAT_6mA_dorado_5khz_v4sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_sup_5kHz_v4_6mA.bed",sep = '\t', header = T)
HP26695_WGA_6mA_dorado_5khz_v4sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_sup_5kHz_v4_6mA.bed",sep = '\t', header = T)

HP26695_NAT_6mA_dorado_5khz_v4hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_hac_5kHz_v4_6mA.bed",sep = '\t', header = T)
HP26695_WGA_6mA_dorado_5khz_v4hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_hac_5kHz_v4_6mA.bed",sep = '\t', header = T)

HP26695_NAT_6mA_dorado_5khz_v5sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_sup_5kHz_v5_6mA.bed",sep = '\t', header = T)
HP26695_WGA_6mA_dorado_5khz_v5sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_sup_5kHz_v5_6mA.bed",sep = '\t', header = T)

HP26695_NAT_6mA_dorado_5khz_v5hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_hac_5kHz_v5_6mA.bed",sep = '\t', header = T)
HP26695_WGA_6mA_dorado_5khz_v5hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_hac_5kHz_v5_6mA.bed",sep = '\t', header = T)


HP26695_NAT_6mA_dorado_4khz_r10 <- HP26695_NAT_6mA_dorado_4khz_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_6mA_dorado_4khz_r10 <- HP26695_WGA_6mA_dorado_4khz_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_6mA_dorado_5khz_v4sup_r10 <- HP26695_NAT_6mA_dorado_5khz_v4sup_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_6mA_dorado_5khz_v4sup_r10 <- HP26695_WGA_6mA_dorado_5khz_v4sup_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_6mA_dorado_5khz_v4hac_r10 <- HP26695_NAT_6mA_dorado_5khz_v4hac_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_6mA_dorado_5khz_v4hac_r10 <- HP26695_WGA_6mA_dorado_5khz_v4hac_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_6mA_dorado_5khz_v5sup_r10 <- HP26695_NAT_6mA_dorado_5khz_v5sup_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_6mA_dorado_5khz_v5sup_r10 <- HP26695_WGA_6mA_dorado_5khz_v5sup_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_6mA_dorado_5khz_v5hac_r10 <- HP26695_NAT_6mA_dorado_5khz_v5hac_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_6mA_dorado_5khz_v5hac_r10 <- HP26695_WGA_6mA_dorado_5khz_v5hac_r10 %>% mutate(sample="HP26695_WGA")


HP26695_6mA_r10_meta <- rbind(HP26695_NAT_6mA_dorado_4khz_r10,HP26695_WGA_6mA_dorado_4khz_r10,
                              HP26695_NAT_6mA_dorado_5khz_v4sup_r10,HP26695_WGA_6mA_dorado_5khz_v4sup_r10,
                              HP26695_NAT_6mA_dorado_5khz_v4hac_r10,HP26695_WGA_6mA_dorado_5khz_v4hac_r10,
                              HP26695_NAT_6mA_dorado_5khz_v5sup_r10,HP26695_WGA_6mA_dorado_5khz_v5sup_r10,
                              HP26695_NAT_6mA_dorado_5khz_v5hac_r10,HP26695_WGA_6mA_dorado_5khz_v5hac_r10)

################ Dorado 4mC dataframes #################

HP26695_NAT_4mC_dorado_5khz_v4sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_sup_5kHz_v4_4mC.bed",sep = '\t', header = T)
HP26695_WGA_4mC_dorado_5khz_v4sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_sup_5kHz_v4_4mC.bed",sep = '\t', header = T)

HP26695_NAT_4mC_dorado_5khz_v5sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_sup_5kHz_v5_4mC.bed",sep = '\t', header = T)
HP26695_WGA_4mC_dorado_5khz_v5sup_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_sup_5kHz_v5_4mC.bed",sep = '\t', header = T)

HP26695_NAT_4mC_dorado_5khz_v5hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_NAT_rep01_kit14_5kHz_hac_5kHz_v5_4mC.bed",sep = '\t', header = T)
HP26695_WGA_4mC_dorado_5khz_v5hac_r10 <- read.csv(file = "HP26695_meta_20x/HP26695_WT_WGA_rep01_kit14_5kHz_hac_5kHz_v5_4mC.bed",sep = '\t', header = T)


HP26695_NAT_4mC_dorado_5khz_v4sup_r10 <- HP26695_NAT_4mC_dorado_5khz_v4sup_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_4mC_dorado_5khz_v4sup_r10 <- HP26695_WGA_4mC_dorado_5khz_v4sup_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_4mC_dorado_5khz_v5sup_r10 <- HP26695_NAT_4mC_dorado_5khz_v5sup_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_4mC_dorado_5khz_v5sup_r10 <- HP26695_WGA_4mC_dorado_5khz_v5sup_r10 %>% mutate(sample="HP26695_WGA")

HP26695_NAT_4mC_dorado_5khz_v5hac_r10 <- HP26695_NAT_4mC_dorado_5khz_v5hac_r10 %>% mutate(sample="HP26695_NAT")
HP26695_WGA_4mC_dorado_5khz_v5hac_r10 <- HP26695_WGA_4mC_dorado_5khz_v5hac_r10 %>% mutate(sample="HP26695_WGA")


HP26695_4mC_r10_meta <- rbind(HP26695_NAT_4mC_dorado_5khz_v4sup_r10,HP26695_WGA_4mC_dorado_5khz_v4sup_r10,
                              HP26695_NAT_4mC_dorado_5khz_v5sup_r10,HP26695_WGA_4mC_dorado_5khz_v5sup_r10,
                              HP26695_NAT_4mC_dorado_5khz_v5hac_r10,HP26695_WGA_4mC_dorado_5khz_v5hac_r10)


################ Deepmod 5mC dataframes #################

HP26695_NAT_deepmod_4khz_r10_bilstm <- read.csv(file = "HP26695_WT_NAT_rep01_kit14_4kHz_100x_bilstmv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
HP26695_NAT_deepmod_5khz_r10_bilstm <- read.csv(file = "HP26695_WT_NAT_rep01_kit14_5kHz_100x_bilstmv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
HP26695_WGA_deepmod_4khz_r10_bilstm <- read.csv(file = "HP26695_WT_WGA_rep01_kit14_4kHz_100x_bilstmv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
HP26695_WGA_deepmod_5khz_r10_bilstm <- read.csv(file = "HP26695_WT_WGA_rep01_kit14_5kHz_100x_bilstmv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)

HP26695_NAT_deepmod_4khz_r10_transformer <- read.csv(file = "HP26695_WT_NAT_rep01_kit14_4kHz_100x_transformerv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
HP26695_NAT_deepmod_5khz_r10_transformer <- read.csv(file = "HP26695_WT_NAT_rep01_kit14_5kHz_100x_transformerv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
HP26695_WGA_deepmod_4khz_r10_transformer <- read.csv(file = "HP26695_WT_WGA_rep01_kit14_4kHz_100x_transformerv41_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)
HP26695_WGA_deepmod_5khz_r10_transformer <- read.csv(file = "HP26695_WT_WGA_rep01_kit14_5kHz_100x_transformerv43_deepmod2_out/output.per_site.withFullContext.tsv",sep = '\t', header = T)


HP26695_NAT_deepmod_4khz_r10_bilstm <- HP26695_NAT_deepmod_4khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="4kHz",sample="HP26695_NAT",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_NAT_deepmod_5khz_r10_bilstm <- HP26695_NAT_deepmod_5khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="5kHz",sample="HP26695_NAT",device="grid",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_WGA_deepmod_4khz_r10_bilstm <- HP26695_WGA_deepmod_4khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="4kHz",sample="HP26695_WGA",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_WGA_deepmod_5khz_r10_bilstm <- HP26695_WGA_deepmod_5khz_r10_bilstm %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="bilstm",sample_rate="5kHz",sample="HP26695_WGA",device="grid",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_NAT_deepmod_4khz_r10_transformer <- HP26695_NAT_deepmod_4khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="4kHz",sample="HP26695_NAT",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_NAT_deepmod_5khz_r10_transformer <- HP26695_NAT_deepmod_5khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="5kHz",sample="HP26695_NAT",device="grid",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_WGA_deepmod_4khz_r10_transformer <- HP26695_WGA_deepmod_4khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="4kHz",sample="HP26695_WGA",device="prom",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_WGA_deepmod_5khz_r10_transformer <- HP26695_WGA_deepmod_5khz_r10_transformer %>% rename("chrom"=X.chromosome,"p1"=position_before,"p2"=position,"per"=mod_fraction,"M"=mod_coverage,"UM"=unmod_coverage) %>%
  mutate(model="transformer",sample_rate="5kHz",sample="HP26695_WGA",device="grid",mod="5mC",flowcell="r10.4.1",tool="deepmod2",rep="rep01",species="HP26695") %>%
  select("chrom","p1","p2","mod","coverage","strand","M","UM","device","flowcell","tool","model","sample","rep","sample_rate","full_context","species","per")

HP26695_deepmod_r10 <- rbind(HP26695_NAT_deepmod_4khz_r10_bilstm,HP26695_WGA_deepmod_4khz_r10_bilstm,
                            HP26695_NAT_deepmod_5khz_r10_bilstm,HP26695_WGA_deepmod_5khz_r10_bilstm,
                            HP26695_NAT_deepmod_4khz_r10_transformer,HP26695_WGA_deepmod_4khz_r10_transformer,
                            HP26695_NAT_deepmod_5khz_r10_transformer,HP26695_WGA_deepmod_5khz_r10_transformer)


HP26695_allTools <- rbind(HP26695_5mC_r10_meta,HP26695_6mA_r10_meta,HP26695_4mC_r10_meta,HP26695_deepmod_r10)

write.table(HP26695_allTools, file = "HP26695_allTools_meta.tsv",sep ="\t",row.names=FALSE,quote = F)


