library(tidyverse)
library(cowplot)
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
  axis.text.x = element_text(angle = 45, hjust = 1),
  strip.text = element_text(size = 24)
)

colors <- c("Bisulfite" = "#66C5CC",
            "DeepMod2_4kHz_BiLSTM" = "#F6CF71",
            "Dorado_4kHz_v4" = "#F89C74",
            "DeepMod2_5kHz_BiLSTM" = "#DCB0F2",
            "Dorado_5kHz_v4" = "#9EB9F3",
            "Dorado_5kHz_v5" = "#FE88B1")


########################################################
# Ecoli AXCG Plot
########################################################

Ecoli_allTools <- read.csv(file = "Ecoli_allTools_meta.tsv",sep = '\t', header = T)

Ecoli_6mA_r10_ACG <- Ecoli_allTools  %>%  filter(grepl(".....ACG.........", full_context),mod=="6mA") %>%
  mutate(context="ACG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% 
  filter(coverage>20) %>% mutate(subcontext=str_sub(full_context,6, 7)) 

Ecoli_6mA_r10_ANCG <- Ecoli_allTools  %>%  filter(grepl(".....A.CG........", full_context),mod=="6mA") %>%
  mutate(context="ANCG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% 
  filter(coverage>20) %>% mutate(subcontext=str_sub(full_context,6, 8))

Ecoli_6mA_r10_ANNCG <- Ecoli_allTools  %>%  filter(grepl(".....A..CG.......", full_context),mod=="6mA") %>%
  mutate(context="ANNCG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% 
  filter(coverage>20) %>% mutate(subcontext=str_sub(full_context,6, 9)) %>% filter(!grepl("CG",subcontext))

Ecoli_6mA_r10_ANNNCG <- Ecoli_allTools  %>%  filter(grepl(".....A...CG......", full_context),mod=="6mA") %>%
  mutate(context="ANNNCG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% 
  filter(coverage>20) %>% mutate(subcontext=str_sub(full_context,6, 10)) %>% filter(!grepl("CG",subcontext))

Ecoli_6mA_r10_ANNNNCG <- Ecoli_allTools  %>%  filter(grepl(".....A....CG.....", full_context),mod=="6mA") %>%
  mutate(context="ANNNNCG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% 
  filter(coverage>20) %>% mutate(subcontext=str_sub(full_context,6, 11)) %>% filter(!grepl("CG",subcontext))

Ecoli_6mA_r10_ANNNNNCG <- Ecoli_allTools  %>%  filter(grepl(".....A.....CG....", full_context),mod=="6mA") %>%
  mutate(context="ANNNNNCG",legend = if_else(tool != "bisulfite",paste(tool,sample_rate,model,sep="_"),paste(tool,sep="_"))) %>% 
  filter(coverage>20) %>% mutate(subcontext=str_sub(full_context,6, 12)) %>% filter(!grepl("CG",subcontext))

Ecoli_6mA_r10_AXCG <- rbind(Ecoli_6mA_r10_ACG,Ecoli_6mA_r10_ANCG,Ecoli_6mA_r10_ANNCG,Ecoli_6mA_r10_ANNNCG,Ecoli_6mA_r10_ANNNNCG,Ecoli_6mA_r10_ANNNNNCG) %>% 
  mutate(legend = case_when(legend=="bisulfite"~"Bisulfite",legend=="deepmod2_4kHz_bilstm"~"DeepMod2_4kHz_BiLSTM",
                            legend=="deepmod2_5kHz_bilstm"~"DeepMod2_5kHz_BiLSTM",legend=="dorado_v4_4kHz_v4_sup"~"Dorado_4kHz_v4",
                            legend=="dorado_v4_5kHz_v4_sup"~"Dorado_5kHz_v4",legend=="dorado_v5_5kHz_v5_sup"~"Dorado_5kHz_v5",TRUE~legend))


Ecoli_6mA_r10_AXCG$legend <- factor(Ecoli_6mA_r10_AXCG$legend,levels = c("Dorado_4kHz_v4","Dorado_5kHz_v4","Dorado_5kHz_v5"))
Ecoli_6mA_r10_AXCG <- Ecoli_6mA_r10_AXCG %>% drop_na()

Ecoli_6mA_r10_AXCG  %>% filter(sample==c("DM","DM+MSssI")) %>%
  ggplot(aes(factor(context), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  facet_wrap(vars(sample))+
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli", fill = "Model",color= "Model")

showtext_opts(dpi = 300)
ggsave('Ecoli_AXCG.png', dpi = 300, bg = "white", height = 6, width = 12)  
showtext_opts(dpi = 96)
showtext_auto()



Ecoli_6mA_r10_AXCG  %>% filter(sample==c("DM+MSssI")) %>%
  ggplot(aes(factor(context), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  facet_wrap(vars(sample))+
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_MSssI_AXCG.png', dpi = 300, bg = "white", height = 6, width = 12)  
showtext_opts(dpi = 96)
showtext_auto()


Ecoli_6mA_r10_AXCG  %>% filter(sample==c("DM")) %>%
  ggplot(aes(factor(context), per))+
  geom_half_boxplot(aes(fill=legend),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) +
  geom_half_point(aes(color=legend),side = "l",shape = ".") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  facet_wrap(vars(sample))+
  ylab("% Methylation")+
  xlab("") +
  theme_minimal()+
  my_theme+
  labs(title = "Species: Ecoli", fill = "Model",color= "Model")


showtext_opts(dpi = 300)
ggsave('Ecoli_DM_AXCG.png', dpi = 300, bg = "white", height = 6, width = 12)  
showtext_opts(dpi = 96)
showtext_auto()


