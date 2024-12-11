########################
# Reading the input file
########################

data<-"performance_metrics.csv"
df<-read.csv(data, header=T)
View(df)

# The header has the following columns
# True negative, False positive, True positive, False negative, Precision, Recall, F1, Specificity, Cutoff, Model, Method, Depth_cutoff, Locations, Species, Motif, Mod, Tool


#################################
# Loading the required libraries
#################################

library(ggplot2)
library(tidyverse)
library(reshape2)
library(showtext)
library(ggtext)

# Setting the theme
my_theme <- theme(
  text = element_text(family = "lato", size = 14),
  plot.title = element_text(size = 14),
  axis.text = element_text(size = 14),
  axis.title = element_text(size=16),
  axis.ticks= element_blank(),
  legend.text = element_text(size=14),
  legend.title = element_text(size=16),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.text = element_text(size = 14),
  panel.grid.major = element_line(colour="gray95")
)


################################################
# Plotting the performance metrics as bar plots
################################################

array <- c("F1", "Precision", "Recall", "Specificity")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting the non-CG 5mC E.coli data (CmCWGG and GmC motifs) as a bar plot

for (metric in array) {

	plot_non_CG <- df %>% filter(Motif %in% c("GmC","CmCWGG") & Cutoff %in% "95-0") %>% ggplot(aes(x=Motif, y=.data[[metric]], fill=Model)) + geom_bar(stat="identity", width=0.7, position=position_dodge(0.8)) + labs(title=paste("Comparison of ",metric," values (non-CG motifs)| E.coli", sep=""), x="", y=paste(metric), fill="Model") + scale_fill_manual(values=c("#F89C74", "#9EB9F3","#FE88B1"),labels=c("Dorado_4kHz_v4","Dorado_5kHz_v4","Dorado_5kHz_v5")) + scale_x_discrete(limits=c("CmCWGG","GmC")) + my_theme + theme(axis.text.x = element_text(colour="black"))
	plot_non_CG
	
	# Saving the plot as a png file
	showtext_opts(dpi = 300)
	ggsave(paste("E.coli_non-CG_",metric,"_barplot.png", sep=""), dpi = 300, bg = "white", height = 4, width = 6)  
	showtext_opts(dpi = 96)
	showtext_auto()
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting the CG 5mC E.coli data as a bar plot

df_CG<-subset(df, Motif=="mCG" & Cutoff=="100-0")
View(df_CG)


df_melt_CG = melt(df_CG, id.vars = c('True.negative','False.positive','True.positive','False.negative','Cutoff','Model','Method','Depth_cutoff','Locations','Species','Motif','Mod','Tool','Rebase'))
View(df_melt_CG)

plot_CG<-ggplot(df_melt_CG, aes(x=Motif, y=value, fill=Model)) + geom_bar(stat="identity", width=0.8, position=position_dodge(0.9)) + labs(title="Comparison of F1 scores (CG motif)| E.coli", x="", y="Score", fill="Model") + scale_fill_manual(values=c("#F6CF71","#F89C74","#DCB0F2","#9EB9F3","#FE88B1"),labels=c("DeepMod2_4kHz","Dorado_4kHz_v4","DeepMod2_5kHz","Dorado_5kHz_v4","Dorado_5kHz_v5")) + my_theme + facet_wrap(~ factor(variable, levels=c("F1","Precision","Recall","Specificity")), strip.position="bottom", nrow=1) + theme(strip.background=element_blank(), axis.text.x=element_blank())
plot_CG

# Saving the plot as a png file
showtext_opts(dpi = 300)
ggsave('E.coli_CG_barplot.png', dpi = 300, bg = "white", height = 3.5, width = 8)  
showtext_opts(dpi = 96)
showtext_auto()


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting the E.coli and H.pylori 6mA data for the various motifs as a bar plot

df_6mA<-subset(df, Mod=="6mA")
View(df_6mA)

my_theme_6mA <- theme(
	text = element_text(family = "lato", size = 14), 
	plot.title = element_text(size = 14), 
	axis.text = element_text(size = 13.5), 
	axis.title = element_text(size=16), 
	axis.ticks= element_blank(), 
	legend.text = element_text(size=12), 
	legend.title = element_text(size=14), 
	panel.grid.minor = element_blank(), 
	panel.background = element_blank(), 
	strip.text = element_text(size = 14), 
	panel.grid.major = element_line(colour="gray95"), 
	axis.text.x=element_text(hjust=0.9, angle=45,size=10)
)


for (metric in array) {

	plot_6mA <-ggplot(df_6mA, aes(x=Motif, y=.data[[metric]], fill=Model)) + geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) + labs(title=paste("Comparison of ",metric," values | 6mA",sep=""), x="", y=paste(metric), fill="Model") + scale_fill_manual(values=c("#F89C74", "#9EB9F3","#FE88B1"),labels=c("Dorado_4kHz_v4","Dorado_5kHz_v4","Dorado_5kHz_v5")) + my_theme_6mA + facet_grid(~ Species, scales="free_x", space="free_x") + theme(strip.background = element_blank(), panel.spacing=unit(1,"lines"))
	plot_6mA
	
	# Saving the plot as a png file
	showtext_opts(dpi = 300)
	ggsave(paste("6mA_",metric,"_barplot.png", sep=""), dpi = 300, bg = "white", height = 4, width = 6.5)  
	showtext_opts(dpi = 96)
	showtext_auto()
}


