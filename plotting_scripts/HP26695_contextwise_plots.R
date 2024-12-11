library(tidyverse)
library(ggplot2)
library(gghalves)
library(Hmisc)
library(paletteer)
library(corrr)
library(reshape2)
library(showtext)


my_theme <- theme(
  text = element_text(family = "roboto", size = 20),
  axis.text = element_text(size = 24),
  axis.title = element_text(size=24),
  legend.text = element_text(size=22),
  legend.title = element_text(size=24),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  legend.key.spacing.y = unit(0.3, "cm"),
  legend.key.spacing = unit(0.3, "cm"),
  strip.text = element_text(size = 24)
)

colors <- c("Bisulfite" = "#66C5CC",
            "DeepMod2_4kHz_BiLSTM" = "#F6CF71",
            "Dorado_4kHz_v4" = "#F89C74",
            "DeepMod2_5kHz_BiLSTM" = "#DCB0F2",
            "Dorado_5kHz_v4" = "#9EB9F3",
            "Dorado_5kHz_v5" = "#FE88B1")


HP26695_allTools <- read.csv(file = "HP26695_allTools_meta.tsv",sep = '\t', header = T)

########################################################
# Generating 6mA context plots 
########################################################

################ GATC context plot #################

HP26695_allTools_GATC <- HP26695_allTools  %>%  filter(grepl("....GATC.........", full_context), mod=="6mA") %>% 
  mutate(context="GATC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_GATC <- HP26695_allTools_GATC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_GATC$legend <- factor(HP26695_allTools_GATC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_GATC$sample <- factor(HP26695_allTools_GATC$sample,levels = c("WGA","NAT"))
HP26695_allTools_GATC <- HP26695_allTools_GATC %>% drop_na()


HP26695_allTools_GATC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: GaTC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_GaTC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ CATG context plot #################

HP26695_allTools_CATG <- HP26695_allTools  %>%  filter(grepl("....CATG.........", full_context), mod=="6mA") %>% 
  mutate(context="CATG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_CATG <- HP26695_allTools_CATG %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_CATG$legend <- factor(HP26695_allTools_CATG$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_CATG$sample <- factor(HP26695_allTools_CATG$sample,levels = c("WGA","NAT"))
HP26695_allTools_CATG <- HP26695_allTools_CATG %>% drop_na()


HP26695_allTools_CATG %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: CaTG", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_CaTG.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ ATTAAT context plot #################

HP26695_allTools_ATTAAT <- HP26695_allTools  %>%  filter(grepl(".ATTAAT..........", full_context), mod=="6mA") %>% 
  mutate(context="ATTAAT",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_ATTAAT <- HP26695_allTools_ATTAAT %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_ATTAAT$legend <- factor(HP26695_allTools_ATTAAT$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_ATTAAT$sample <- factor(HP26695_allTools_ATTAAT$sample,levels = c("WGA","NAT"))
HP26695_allTools_ATTAAT <- HP26695_allTools_ATTAAT %>% drop_na()


HP26695_allTools_ATTAAT %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: ATTAaT", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_ATTAaT.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ TCGA context plot #################

HP26695_allTools_TCGA <- HP26695_allTools  %>%  filter(grepl("..TCGA...........", full_context), mod=="6mA") %>% 
  mutate(context="TCGA",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_TCGA <- HP26695_allTools_TCGA %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_TCGA$legend <- factor(HP26695_allTools_TCGA$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_TCGA$sample <- factor(HP26695_allTools_TCGA$sample,levels = c("WGA","NAT"))
HP26695_allTools_TCGA <- HP26695_allTools_TCGA %>% drop_na()


HP26695_allTools_TCGA %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: TCGa", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_TCGa.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ GCAG context plot #################

HP26695_allTools_GCAG <- HP26695_allTools  %>%  filter(grepl("...GCAG..........", full_context), mod=="6mA") %>% 
  mutate(context="GCAG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_GCAG <- HP26695_allTools_GCAG %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_GCAG$legend <- factor(HP26695_allTools_GCAG$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_GCAG$sample <- factor(HP26695_allTools_GCAG$sample,levels = c("WGA","NAT"))
HP26695_allTools_GCAG <- HP26695_allTools_GCAG %>% drop_na()


HP26695_allTools_GCAG %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: GCaG", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_GCaG.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ GAAGA context plot #################

HP26695_allTools_GAAGA <- HP26695_allTools  %>%  filter(grepl(".GAAGA...........", full_context), mod=="6mA") %>% 
  mutate(context="GAAGA",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_GAAGA <- HP26695_allTools_GAAGA %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_GAAGA$legend <- factor(HP26695_allTools_GAAGA$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_GAAGA$sample <- factor(HP26695_allTools_GAAGA$sample,levels = c("WGA","NAT"))
HP26695_allTools_GAAGA <- HP26695_allTools_GAAGA %>% drop_na()


HP26695_allTools_GAAGA %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: GAAGa", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_GAAGa.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ GANTC context plot #################

HP26695_allTools_GANTC <- HP26695_allTools  %>%  filter(grepl("....GA.TC........", full_context), mod=="6mA") %>% 
  mutate(context="GANTC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_GANTC <- HP26695_allTools_GANTC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_GANTC$legend <- factor(HP26695_allTools_GANTC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_GANTC$sample <- factor(HP26695_allTools_GANTC$sample,levels = c("WGA","NAT"))
HP26695_allTools_GANTC <- HP26695_allTools_GANTC %>% drop_na()


HP26695_allTools_GANTC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 6mA |  Context: GaNTC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_GaNTC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()


########################################################
# Generating 5mC context plots 
########################################################

################ CCTC context plot #################

HP26695_allTools_CCTC <- HP26695_allTools  %>%  filter(grepl(".....CCTC........", full_context), mod=="5mC") %>% 
  mutate(context="CCTC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_CCTC <- HP26695_allTools_CCTC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_CCTC$legend <- factor(HP26695_allTools_CCTC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_CCTC$sample <- factor(HP26695_allTools_CCTC$sample,levels = c("WGA","NAT"))
HP26695_allTools_CCTC <- HP26695_allTools_CCTC %>% drop_na()


HP26695_allTools_CCTC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 5mC |  Context: cCTC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_cCTC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ GCGC context plot #################

HP26695_allTools_GCGC <- HP26695_allTools  %>%  filter(grepl("....GCGC.........", full_context), mod=="5mC") %>% 
  mutate(context="GCGC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_GCGC <- HP26695_allTools_GCGC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_GCGC$legend <- factor(HP26695_allTools_GCGC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_GCGC$sample <- factor(HP26695_allTools_GCGC$sample,levels = c("WGA","NAT"))
HP26695_allTools_GCGC <- HP26695_allTools_GCGC %>% drop_na()


HP26695_allTools_GCGC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 5mC |  Context: GcGC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_GcGC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()


########################################################
# Generating 4mC context plots 
########################################################

################ TCTTC context plot #################

HP26695_allTools_TCTTC <- HP26695_allTools  %>%  filter(grepl("....TCTTC........", full_context), mod=="4mC") %>% 
  mutate(context="TCTTC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
HP26695_allTools_TCTTC <- HP26695_allTools_TCTTC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_4kHz_sup_v4"~"Dorado_4kHz_v4",
                            legend=="dorado_5kHz_sup_v4"~"Dorado_5kHz_v4",legend=="dorado_5kHz_sup_v5"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="HP26695_WGA"~"WGA",sample=="HP26695_NAT"~"NAT", TRUE~sample))
HP26695_allTools_TCTTC$legend <- factor(HP26695_allTools_TCTTC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
HP26695_allTools_TCTTC$sample <- factor(HP26695_allTools_TCTTC$sample,levels = c("WGA","NAT"))
HP26695_allTools_TCTTC <- HP26695_allTools_TCTTC %>% drop_na()


HP26695_allTools_TCTTC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: HP26695 | Mod: 4mC |  Context: TcTTC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('HP26695_TcTTC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()