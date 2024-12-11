#########################
# Reading the input file
#########################

data <- "minQual_metric.csv"
df_data<-read.csv(data, header=TRUE)
View(df_data)
# The header has the following columns
# True negative, False positive, True positive, False negative, Precision, Recall, F1, Specificity, minQual, Context, Sampling, rate, Model, Mode

#################################
# Loading the required libraries
#################################
library(ggplot2)
library(tidyverse)
library(ggtext)
library(showtext)
library(reshape2)
library(patchwork)


#########################
# Creating the dataframe 
#########################

df_melt<-melt(df_data, id.vars=c("minQual","Context","Sampling.rate","Model","Mode","Plot","True.negative","True.positive","False.negative","False.positive"))
View(df_melt)

# Setting the breaks for the x axis (minimum read quality)
my_breaks<-seq(8,22,2)

# Setting the theme for plotting
my_theme <- theme(
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
	axis.text.x=element_text(hjust=0.5)
)


######################################################################################
# Generating the line plot for E.coli data based on the minimum read quality cut-offs
######################################################################################
 
# The performance metrics are plotted using different linetypes and coloured by the basecalling model across the different minimum read quality datasets.
 
#-------------------------------------------------------------------------------------------------------------------------------------------
# Plotting the line plot for E.coli 5mC wildtype data
WT<-df_melt %>% filter(Context=="WT") %>% ggplot(aes(x=minQual, y=value, group=interaction(Sampling.rate, Context, Model,variable), linetype=variable)) + geom_line(aes(colour=interaction(Sampling.rate,Model,Mode))) + geom_point(size=0.7, colour="grey33") + scale_x_continuous(breaks=my_breaks) + scale_linetype_manual(name = "", values=c("Precision"="dashed","F1"="solid","Specificity"="dotted","Recall"="twodash"), limits=c("F1","Recall","Precision","Specificity")) + labs(title="CmCWGG", x="Minimum read quality", y="Score", color="Model") + scale_color_manual(limits=c("4kHz.v4.sup","5kHz.v4.sup","5kHz.v5.sup"),labels=c("Dorado_4kHz_v4","Dorado_5kHz_v4","Dorado_5kHz_v5"),values=c("#F89C74","#9EB9F3","#FE88B1")) + my_theme + theme(panel.grid.major = element_line(colour="gray95")) + scale_y_continuous(limits=c(0.95,1)) + guides(size=FALSE)
WT

#-------------------------------------------------------------------------------------------------------------------------------------------
# Plotting the line plot for E.coli 5mC MSssI treated data.
MSssI<-df_melt %>% filter(Context=="MSssI") %>% ggplot(aes(x=minQual, y=value, group=interaction(Sampling.rate, Context, Model,variable), linetype=variable)) + geom_line(aes(colour=interaction(Sampling.rate,Model,Mode))) + geom_point(size=0.7, colour="grey33") + scale_x_continuous(breaks=my_breaks) + scale_color_manual(limits=c("4kHz.deepmod.sup","4kHz.v4.sup","5kHz.deepmod.sup","5kHz.v4.sup","5kHz.v5.sup"),values=c("#F6CF71","#F89C74","#DCB0F2","#9EB9F3","#FE88B1"),labels=c("DeepMod2_4kHz","Dorado_4kHz_v4","DeepMod2_5kHz","Dorado_5kHz_v4","Dorado_5kHz_v5")) + scale_linetype_manual(name = "", values=c("Precision"="dashed","F1"="solid","Specificity"="dotted","Recall"="twodash"), limits=c("F1","Recall","Precision","Specificity")) + labs(title="mCG", x="Minimum read quality", y="Score", color="Model") + my_theme + theme(panel.grid.major = element_line(colour="gray95")) + scale_y_continuous(limits=c(0.90,1)) + guides(size=FALSE) + theme(axis.text = element_text(size = 14), plot.title = element_text(size = 14.5))
MSssI

#-------------------------------------------------------------------------------------------------------------------------------------------
# Plotting the line plot for E.coli 6mA wildtype data.
GATC<-df_melt %>% filter(Context=="6mA") %>% ggplot(aes(x=minQual, y=value, group=interaction(Sampling.rate, Context, Model,variable), linetype=variable)) + geom_line(aes(colour=interaction(Sampling.rate,Model,Mode))) + geom_point(size=0.7, colour="grey33") + scale_x_continuous(breaks=my_breaks) + scale_linetype_manual(name = "", values=c("Precision"="dashed","F1"="solid","Specificity"="dotted","Recall"="twodash"), limits=c("F1","Recall","Precision","Specificity")) + labs(title="GmATC", x="Minimum read quality", y="Score", color="Model") + scale_color_manual(limits=c("4kHz.v4.sup","5kHz.v4.sup","5kHz.v5.sup"),labels=c("Dorado_4kHz_v4","Dorado_5kHz_v4","Dorado_5kHz_v5"),values=c("#F89C74","#9EB9F3","#FE88B1")) + my_theme + theme(panel.grid.major = element_line(colour="gray95")) + scale_y_continuous(limits=c(0.95,1)) + guides(size=FALSE)
GATC

#-------------------------------------------------------------------------------------------------------------------------------------------

##########################################
# Combining all the individual line plots
##########################################

((GATC + WT + plot_layout(guides="collect", axis_titles = "collect")) / ((MSssI + plot_spacer()) + plot_layout(widths=c(2,0.5))) ) + plot_annotation('Effect of minimum read quality', theme=theme(plot.title = element_text(size=15)))

# Saving the line plots as a png file
showtext_opts(dpi = 300)
ggsave('minQual_effect_patchwork.png', dpi = 300, bg = "white", height = 7, width = 12)  
showtext_opts(dpi = 96)
showtext_auto()

