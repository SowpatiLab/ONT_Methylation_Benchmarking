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


Ecoli_allTools <- read.csv(file = "Ecoli_allTools_meta.tsv",sep = '\t', header = T)

########################################################
# Generating 5mC context plots 
########################################################

################ CG context plot #################

# Samples : "Ecoli_DM","Ecoli_MSssI"
Ecoli_allTools_CG <- Ecoli_allTools  %>%  filter(grepl(".....CG..........", full_context),mod=="5mC") %>%
  mutate(context="CG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
Ecoli_allTools_CG <- Ecoli_allTools_CG %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_v4_4kHz_v4_sup"~"Dorado_4kHz_v4",
                            legend=="dorado_v4_5kHz_v4_sup"~"Dorado_5kHz_v4",legend=="dorado_v5_5kHz_v5_sup"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="Ecoli_DM"~"DM",sample=="Ecoli_MSssI"~"MSssI", TRUE~sample))
Ecoli_allTools_CG$legend <- factor(Ecoli_allTools_CG$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
Ecoli_allTools_CG$sample <- factor(Ecoli_allTools_CG$sample,levels = c("DM","MSssI"))
Ecoli_allTools_CG <- Ecoli_allTools_CG %>% drop_na()


Ecoli_allTools_CG %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli | Mod: 5mC |  Context: cG", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_cG.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ CCWGG context plot #################

# Samples : "Ecoli_DM","Ecoli_WT"
Ecoli_allTools_CCWGG <- Ecoli_allTools  %>%  filter(grepl("....CC[AT]GG........", full_context),mod=="5mC") %>%
  mutate(context="CCWGG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
Ecoli_allTools_CCWGG <- Ecoli_allTools_CCWGG %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                         legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_v4_4kHz_v4_sup"~"Dorado_4kHz_v4",
                           legend=="dorado_v4_5kHz_v4_sup"~"Dorado_5kHz_v4",legend=="dorado_v5_5kHz_v5_sup"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="Ecoli_DM"~"DM",sample=="Ecoli_WT"~"WT", TRUE~sample))

Ecoli_allTools_CCWGG$legend <- factor(Ecoli_allTools_CCWGG$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
Ecoli_allTools_CCWGG$sample <- factor(Ecoli_allTools_CCWGG$sample,levels = c("DM","WT"))
Ecoli_allTools_CCWGG <- Ecoli_allTools_CCWGG %>% drop_na()

Ecoli_allTools_CCWGG %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli | Mod: 5mC |  Context: CcWGG", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_CcWGG.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ GC context plot #################

# Samples : "Ecoli_DM","Ecoli_MCviP"
Ecoli_allTools_GC <- Ecoli_allTools  %>%  filter(grepl("....GC...........", full_context),mod=="5mC") %>%
  mutate(context="GC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
Ecoli_allTools_GC <- Ecoli_allTools_GC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_v4_4kHz_v4_sup"~"Dorado_4kHz_v4",
                            legend=="dorado_v4_5kHz_v4_sup"~"Dorado_5kHz_v4",legend=="dorado_v5_5kHz_v5_sup"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="Ecoli_DM"~"DM",sample=="Ecoli_MCviP"~"MCviP", TRUE~sample))
Ecoli_allTools_GC$legend <- factor(Ecoli_allTools_GC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
Ecoli_allTools_GC$sample <- factor(Ecoli_allTools_GC$sample,levels = c("DM","MCviP"))
Ecoli_allTools_GC <- Ecoli_allTools_GC %>% drop_na()


Ecoli_allTools_GC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli | Mod: 5mC |  Context: Gc", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_Gc.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

########################################################
# Generating 6mA context plots 
########################################################

################ GATC context plot #################

# Samples : "Ecoli_DM","Ecoli_WT"
Ecoli_allTools_GATC <- Ecoli_allTools  %>%  filter(grepl("....GATC.........", full_context),mod=="6mA") %>%
  mutate(context="GATC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
Ecoli_allTools_GATC <- Ecoli_allTools_GATC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_v4_4kHz_v4_sup"~"Dorado_4kHz_v4",
                            legend=="dorado_v4_5kHz_v4_sup"~"Dorado_5kHz_v4",legend=="dorado_v5_5kHz_v5_sup"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="Ecoli_DM"~"DM",sample=="Ecoli_WT"~"WT", TRUE~sample))
Ecoli_allTools_GATC$legend <- factor(Ecoli_allTools_GATC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
Ecoli_allTools_GATC$sample <- factor(Ecoli_allTools_GATC$sample,levels = c("DM","WT"))
Ecoli_allTools_GATC <- Ecoli_allTools_GATC %>% drop_na()


Ecoli_allTools_GATC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli | Mod: 6mA |  Context: GaTC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_GaTC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

################ AACNNNNNNGTGC context plot #################

# Samples : "Ecoli_DM","Ecoli_WT"
Ecoli_allTools_AACNNNNNNGTGC <- Ecoli_allTools  %>%  filter(grepl("....AAC......GTGC", full_context),mod=="6mA") %>%
  mutate(context="AACNNNNNNGTGC",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% filter(coverage>20)
Ecoli_allTools_AACNNNNNNGTGC <- Ecoli_allTools_AACNNNNNNGTGC %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_v4_4kHz_v4_sup"~"Dorado_4kHz_v4",
                            legend=="dorado_v4_5kHz_v4_sup"~"Dorado_5kHz_v4",legend=="dorado_v5_5kHz_v5_sup"~"Dorado_5kHz_v5",TRUE~legend),
         sample= case_when(sample=="Ecoli_DM"~"DM",sample=="Ecoli_WT"~"WT", TRUE~sample))
Ecoli_allTools_AACNNNNNNGTGC$legend <- factor(Ecoli_allTools_AACNNNNNNGTGC$legend,levels = c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"))
Ecoli_allTools_AACNNNNNNGTGC$sample <- factor(Ecoli_allTools_AACNNNNNNGTGC$sample,levels = c("DM","WT"))
Ecoli_allTools_AACNNNNNNGTGC <- Ecoli_allTools_AACNNNNNNGTGC %>% drop_na()


Ecoli_allTools_AACNNNNNNGTGC %>%
  ggplot(aes(factor(sample), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli | Mod: 6mA |  Context: AaCNNNNNNGTGC", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_AaCNNNNNNGTGC.png', dpi = 300, bg = "white", height =6, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()