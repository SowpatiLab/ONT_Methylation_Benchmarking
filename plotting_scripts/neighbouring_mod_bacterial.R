library(dplyr)
library(ggplot2)
library(gghalves)
library(ggtext)
library(patchwork)
library(paletteer)
library(forcats)
library(showtext)
library(stringr)
library(arrow)
showtext_auto()


# Set theme
font_add_google("Roboto", "roboto")
font_add_google("Lato", "lato")
font_add_google("Ubuntu", "ubuntu")
font_add_google("Oxygen", "oxygen")
font_add_google("Noto Sans", 'noto')
font_add_google("Nunito", 'nunito')

my_theme <- theme(
  text = element_text(family = "lato", size = 14),
  plot.title = element_markdown(size = 14),
  axis.text = element_text(size = 14),
  axis.title = element_text(size=16),
  legend.text = element_text(size=14),
  legend.title = element_text(size=16),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.text = element_text(size = 14)
)

# User defined functions
pattern <- function (x, base = 'a', lpad = 5, rpad = 11) {
  pos = str_locate(x, base)[1]
  len = str_length(x)
  prefix = str_sub(x, start = 1, end = pos-1)
  suffix = str_sub(x, start = pos+1, end = len)
  left = lpad - str_length(prefix)
  right = rpad - str_length(suffix)
  final = str_glue(str_dup('.', left), prefix, str_to_upper(base), suffix, str_dup('.', right), sep = '')
  final
}
pattern('ATTAaT')
pattern('ATTAcT', base = 'c')

ignoreG <- function(x, lpad = 5, rpad = 0) {
  start = str_locate_all(x, "\\.+")[[1]][1]
  end = str_locate_all(x, "\\.+")[[1]][2]
  c(start+lpad, end+lpad)
}


colors <- c("Bisulfite" = "#66C5CC",
            "DeepMod2_4kHz_BiLSTM" = "#F6CF71",
            "Dorado_4kHz_v4" = "#F89C74",
            "DeepMod2_5kHz_BiLSTM" = "#DCB0F2",
            "Dorado_5kHz_v4" = "#9EB9F3",
            "Dorado_5kHz_v5" = "#FE88B1")

# Reading data
df_Ecoli <- read.csv("./Ecoli_allTools_meta.tsv", header = T, sep = "\t", stringsAsFactors = F)
head(df_Ecoli)

df_Hpylori <- read.csv("./HP26695_allTools_meta.tsv", header = T, sep = "\t", stringsAsFactors = F)
head(df_Hpylori)


################################################################################
# Species: E coli
# Cleaning up dataframe

df_Ecoli %>%
  count(model, mod)

# Remove sites with <20x coverage
df_Ecoli <- df_Ecoli %>%
  filter(coverage >= 20)

# Remove transfomer models of DeepMod2
df_Ecoli <- df_Ecoli %>%
  filter(model != "transformer")

# Remove high accuracy calls and retain only sup calls
df_Ecoli <- df_Ecoli %>%
  filter(!(grepl("hac", model)))
head(df_Ecoli)

df_Ecoli %>%
  count(model)

# Rename Tool and Condition names for clarity
df_Ecoli <- df_Ecoli %>%
  mutate(condition = case_when(
    model == '-' ~ "Bisulfite",
    model == 'bilstm' ~ str_c("DeepMod2_", sample_rate, "_", "BiLSTM"),
    .default = str_replace(str_c("Dorado_", sample_rate, "_", model), "_sup", "")
  )) %>%
  mutate(name = str_replace(sample, 'Ecoli_', '')) %>%
  mutate(name = case_when(
    name == "MCviP" ~ "DM+M.CviP",
    name == "DM_MCviP" ~ "DM+M.CviP",
    name == "MSssI" ~ "DM+M.SssI",
    name == "DM_MSssI" ~ "DM+M.SssI",
    .default = name
  ))

df_Ecoli %>%
  count(condition, name)

################################################################################
# Species: E coli
# C next to mC

glue_c <- "<span style='color:red'>**C**</span>"
mod <- "5mC"


motif <- "cC.GG"
p1 <- df_Ecoli %>%
# df_Ecoli %>%
  filter(name %in% c("WT", "DM")) %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 8, end = 8) %in% c("A","T"), "CCWGG", "CCSGG")) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("WT", "DM"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("{glue_c}mCNGG"), y = "% Methylation", fill = "Model", title = "Species: *E.coli*")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_{mod}_{str_replace_all(motif, "[.]", "N")}.png'), dpi = 300, bg = "white", height = 3, width = 10)  
showtext_opts(dpi = 96)
showtext_auto() 

s <- "WT"
# s <- "DM"
p1 <- df_Ecoli %>%
  filter(name == s) %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 8, end = 8) %in% c("A","T"), "CCWGG", "CCSGG")) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(~name) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("{glue_c}mCNGG"), y = "% Methylation", fill = "Model", title = str_glue("Species: *E.coli*"))
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_{mod}_{str_replace_all(motif, "[.]", "N")}_{s}.png'), dpi = 300, bg = "white", height = 3, width = 5.5)  
showtext_opts(dpi = 96)
showtext_auto() 


motif <- 'cCG'
p2 <- df_Ecoli %>%
  filter(mod == "5mC") %>%
  filter(name %in% c("DM", "DM+M.SssI")) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  ggplot(aes(x = condition, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(~name) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("{glue_c}mCG"), y = "% Methylation", fill = "Model", title = "Species: *E.coli*") +
  theme(axis.text.x = element_blank())
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_{mod}_{str_replace_all(motif, "[.]", "N")}.png'), dpi = 300, bg = "white", height = 3, width = 5.5)  
showtext_opts(dpi = 96)
showtext_auto() 

p1 + p2 +
  plot_layout(guides = "collect") &
  my_theme
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_Fig4A_combined.png'), dpi = 300, bg = "white", height = 3, width = 8)  
showtext_opts(dpi = 96)
showtext_auto() 


motif <- "CGcG"
df_Ecoli %>%
  filter(name %in% c("DM", "DM+M.SssI")) %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  filter(!grepl('G', str_sub(full_context, start = ignoreG(motif)[1], end = ignoreG(motif)[2]))) %>%
  mutate(seq = str_sub(full_context, start = 6, end = 7)) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("WT", "DM"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("mCG{glue_c}G"), y = "% Methylation", fill = "Model", title = "Species: *E.coli*") +
  theme(axis.text.x = element_blank())
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_{mod}_combined_CGcG.png'), dpi = 300, bg = "white", height = 4, width = 8)  
showtext_opts(dpi = 96)
showtext_auto() 



cutoff <- 50
motif <- "cCG"
motifs <- c("c.CG", "c..CG", "c...CG", "c....CG", "c.....CG", "c......CG", "c.......CG")

motif <- "cCGCG"
motifs <- c("c.CGCG", "c..CGCG", "c...CGCG", "c....CGCG", "c.....CGCG", "c......CGCG", "c.......CGCG")

fp_df <- df_Ecoli %>%
  filter(name == "DM+M.SssI") %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  filter(!grepl('G', str_sub(full_context, start = ignoreG(motif)[1], end = ignoreG(motif)[2]))) %>%
  group_by(condition) %>%
  summarize(total = n(),
            fp = length(condition[per > cutoff]),
            motif = motif,
            distance = 1)

for(motif in motifs) {
  print(motif)
  temp_df <- df_Ecoli %>%
  filter(name == "DM+M.SssI") %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  filter(!grepl('G', str_sub(full_context, start = ignoreG(motif)[1], end = ignoreG(motif)[2]))) %>%
  group_by(condition) %>%
  summarize(total = n(),
            fp = length(condition[per > cutoff]),
            motif = motif,
            distance = 2 + (ignoreG(motif)[2] - ignoreG(motif)[1]))
  fp_df <- rbind(fp_df, temp_df)
}  


p_fpr_cg <- fp_df %>%
  mutate(fpr = fp*100/total) %>%
  ggplot(aes(x = distance, y = fpr, group = condition, colour = condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = "Distance from mCG", y = "False Positive Rate", color = "Model")

p_fpr_cgcg <- fp_df %>%
  mutate(fpr = fp*100/total) %>%
  ggplot(aes(x = distance, y = fpr, group = condition, colour = condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = "Distance from mCGmCG", y = "False Positive Rate", color = "Model")


p_fpr_cg / p_fpr_cgcg +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Species: *E.coli*, Dataset: DM+M.SssI") &
  my_theme
  # theme(legend.position = "bottom")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_5mC_DistancePlot_vertical.png'), dpi = 300, bg = "white", height = 6, width = 5.5)  
showtext_opts(dpi = 96)
showtext_auto() 




################################################################################
# Species: E coli
# A next to mC

df_Ecoli <- df_Ecoli %>%
  mutate(condition = case_when(
    model == '-' ~ "Bisulfite",
    model == 'bilstm' ~ str_c("DeepMod2_", sample_rate, "_", "BiLSTM"),
    .default = str_replace(str_c("Dorado_", sample_rate, "_", model), "_sup", "")
  )) %>%
  mutate(name = sample)
  

glue_a <- "<span style='color:red'>**A**</span>"
mod <- "6mA"

motif <- 'aCG'
p1 <- df_Ecoli %>%
  filter(mod == mod) %>%
  filter(name %in% c("DM", "DM+M.SssI")) %>%
  filter(grepl(pattern(motif, base = 'a'), full_context)) %>%
  ggplot(aes(x = condition, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  # facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  facet_wrap(~name) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("{glue_a}mCG"), y = "% Methylation", fill = "Model", title = "Species: *E.coli*") +
  theme(axis.text.x = element_blank())
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_6mA_aCG.png'), dpi = 300, bg = "white", height = 3, width = 6)  
showtext_opts(dpi = 96)
showtext_auto() 

motif <- 'CCaGG'
p2 <- df_Ecoli %>%
  filter(mod == mod) %>%
  filter(name %in% c("WT", "DM")) %>%
  filter(grepl(pattern(motif, base = 'a'), full_context)) %>%
  ggplot(aes(x = condition, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  # facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  facet_wrap(~fct_relevel(name, c("WT", "DM"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("CmC{glue_a}GG"), y = "% Methylation", fill = "Model", title = "Species: *E.coli*") +
  theme(axis.text.x = element_blank())

showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_6mA_CCaGG.png'), dpi = 300, bg = "white", height = 3, width = 6)  
showtext_opts(dpi = 96)
showtext_auto() 


#Vertical
p2 / p1 +
  plot_layout(guides = "collect") &
  # plot_annotation(title = "Species: *E.coli*, Dataset: DM+M.SssI") &
  my_theme
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_6mA_neighboring_combined.png'), dpi = 300, bg = "white", height = 6, width = 6)  
showtext_opts(dpi = 96)
showtext_auto()

# Horizontal
p2 + p1 +
  plot_layout(guides = "collect") &
  # plot_annotation(title = "Species: *E.coli*, Dataset: DM+M.SssI") &
  my_theme
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Ecoli_6mA_neighboring_combined.png'), dpi = 300, bg = "white", height = 3, width = 8)  
showtext_opts(dpi = 96)
showtext_auto()

  
################################################################################
# Species: H pylori
# Cleaning up the df

df_Hpylori %>%
  count(model)

df_Hpylori <- df_Hpylori %>%
  filter(model != "transformer") %>%
  filter(!(grepl("hac", model)))

df_Hpylori <- df_Hpylori %>%
  mutate(condition = case_when(
    model == '-' ~ "Bisulfite",
    model == 'bilstm' ~ str_c("DeepMod2_", sample_rate, "_", "BiLSTM"),
    .default = str_replace(str_c("Dorado_", sample_rate, "_", model), "_sup", "")
  )) %>%
  mutate(name = str_replace(sample, 'HP26695_', ''))
  

################################################################################
# Species: H pylori
# A next to mC

mod <- "6mA"
glue_a <- "<span style='color:red'>**A**</span>"

motif <- 'aGATC'
df_Hpylori %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'a'), full_context)) %>%
  # mutate(seq = if_else(str_sub(full_context, start = 8, end = 8) %in% c("A","T"), "CCWGG", "CCSGG")) %>%
  ggplot(aes(x = condition, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("{glue_a}GmCGC"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*") +
  theme(axis.text.x = element_blank())



mod <- "5mC"

motif <- "cGCGC"
df_Hpylori %>%
  filter(mod == mod) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  mutate(seq = str_sub(full_context, start = 5, end = 8)) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("GmAG{glue_c}"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*") 
  # theme(axis.text.x = element_blank())


df_Hpylori %>%
  filter(mod == "6mA") %>%
  filter(per < 50) %>%
  ggplot(aes(x = condition, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  labs(x = str_glue("GmAG{glue_c}"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*") 




################################################################################
# Species: H pylori
# mA next to mA

glue_ma <- "<span style='color:red'>**mA**</span>"
glue_mc <- "<span style='color:red'>**mC**</span>"
glue_4mc <- "<span style='color:red'>**4mC**</span>"

s <- "NAT"
# s <- "WGA"
motif <- 'CaTG.TC'
df_Hpylori %>%
  filter(mod == "6mA") %>%
  filter(name == s) %>%
  filter(grepl(pattern(motif, base = 'a'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 9, end = 9) == "A", str_glue("C{glue_ma}TGmATC"), str_glue("C{glue_ma}TGBTC"))) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  theme(axis.text.x = element_markdown()) +
  labs(x = str_glue("C{glue_ma}TGNTC"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Hpylori_{mod}_{str_replace_all(motif, "[.]", "N")}_{s}.png'), dpi = 300, bg = "white", height = 3.5, width = 6)  
showtext_opts(dpi = 96)
showtext_auto()

motif <- 'GaTC.TG'
df_Hpylori %>%
  filter(mod == "6mA") %>%
  filter(name == s) %>%
  filter(grepl(pattern(motif, base = 'a'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 9, end = 9) == "A", str_glue("G{glue_ma}TCmATG"), str_glue("G{glue_ma}TCBTC"))) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  theme(axis.text.x = element_markdown()) +
  labs(x = str_glue("G{glue_ma}TCNTC"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Hpylori_{mod}_{str_replace_all(motif, "[.]", "N")}_{s}.png'), dpi = 300, bg = "white", height = 3.5, width = 6)  
showtext_opts(dpi = 96)
showtext_auto()




motif <- 'GcGCG.'
df_Hpylori %>%
  filter(mod == '5mC') %>%
  filter(name == s) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 10, end = 10) == "C", str_glue("G{glue_mc}GmCGC"), str_glue("G{glue_mc}GCGD"))) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  theme(axis.text.x = element_markdown()) +
  labs(x = str_glue("G{glue_mc}GCGN"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Hpylori_{mod}_{str_replace_all(motif, "[.]", "N")}_{s}.png'), dpi = 300, bg = "white", height = 3.5, width = 7)  
showtext_opts(dpi = 96)
showtext_auto()




motif <- 'GcGC.TC'
df_Hpylori %>%
  filter(mod == '5mC') %>%
  filter(name == s) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 9, end = 9) == "C", str_glue("G{glue_mc}GmCCTC"), str_glue("G{glue_mc}GCDTC"))) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  theme(axis.text.x = element_markdown()) +
  labs(x = str_glue("G{glue_mc}GCNTC"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Hpylori_{mod}_{str_replace_all(motif, "[.]", "N")}_{s}.png'), dpi = 300, bg = "white", height = 3.5, width = 7)  
showtext_opts(dpi = 96)
showtext_auto()



motif <- '.CTcTTC'
df_Hpylori %>%
  filter(mod == '4mC') %>%
  filter(name == s) %>%
  filter(grepl(pattern(motif, base = 'c'), full_context)) %>%
  mutate(seq = if_else(str_sub(full_context, start = 3, end = 3) == "C", str_glue("mCCT{glue_4mc}TTC"), str_glue("DCT{glue_4mc}TTC"))) %>%
  ggplot(aes(x = seq, y = per, fill = condition)) +
  geom_half_boxplot(outlier.shape = NA, lwd = 0.3, fatten = 1) +
  geom_half_point(aes(color = condition), size = 0.1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors, guide = "none") +
  facet_wrap(facets = ~fct_relevel(name, c("NAT", "WGA"))) +
  theme_minimal() +
  my_theme +
  theme(axis.title.x = element_markdown()) +
  theme(axis.text.x = element_markdown()) +
  ylim(c(0,100)) +
  labs(x = str_glue("NCT{glue_4mc}TTC"), y = "% Methylation", fill = "Model", title = "Species: *H.pylori*")
showtext_opts(dpi = 300)
ggsave(str_glue('./plots/neighbors/Hpylori_{mod}_{str_replace_all(motif, "[.]", "N")}_{s}.png'), dpi = 300, bg = "white", height = 3.5, width = 6)  
showtext_opts(dpi = 96)
showtext_auto()








