#################################
# Loading the required libraries
#################################

library(ggplot2)
library(dplyr)
library(gghalves)
library(showtext)
library(ggtext)
library(patchwork)


########################
# Reading the input bed
########################

data_ds_100 <- "df_ds_CpG100.bed"	
# This bed file contains information about the percent ONT and bisulfite methyaltion of cytosines downstream of the fully methylated CpG

data_us_100 <- "df_us_CpG100.bed"
# This bed file contains information about the percent ONT and bisulfite methyaltion of cytosines upstream of the fully methylated CpG

data_ds_0<-"df_ds_CpG0.bed"
# This bed file contains information about the percent ONT and bisulfite methylation of cytosines downstream of the fully unmethylated CpG

data_us_0<-"df_us_CpG0.bed"
# This bed file contains information about the percent ONT and bisulfite methylation of cytosines upstream of the fully unmethylated CpG

ds_100 <- read.csv(data_ds_100, sep="\t", header=F)
us_100 <- read.csv(data_us_100, sep="\t", header=F)

ds_0<-read.csv(data_ds_0, sep="\t", header=F)
us_0<-read.csv(data_us_0, sep="\t", header=F)

# Setting the column names
header <- c("chr", "p1", "p2", "per_ONT", "per_bis", "strand", "index", "model" )

# Setting the theme for plotting
my_theme <- theme(
	text = element_text(family = "lato", size = 14), 
	plot.title = element_text(size = 14), 
	axis.text = element_text(size = 13.5), 
	axis.title = element_text(size=16), 
	axis.ticks= element_blank(), 
	legend.text = element_text(size=14), 
	legend.title = element_text(size=16), 
	panel.grid.minor = element_blank(), 
	panel.background = element_blank(), 
	strip.text = element_text(size = 14), 
	panel.grid.major = element_line(colour="gray95")
)


#####################################################################
# Plotting the neighbouring C methylation status as jitter-box plots
#####################################################################

array <- c("ds_100", "us_100", "ds_0", "us_0")
plot_list <- list()

for (df_name in array) {

	df <- get(df_name)
	names(df) <- header
	View(df)

	# Finding the difference in percent methylation (to plot on the y axis)
	df$difference <- df$per_ONT - df$per_bis

	# Plotting the data 
	plot_list[[df_name]]<-ggplot(df, aes(x=factor(index), y=difference)) + geom_half_boxplot(aes(fill=model),side = "r",errorbar.draw = F,outlier.shape = NA,fatten=1) + geom_half_point(aes(color=model),side = "l",shape = ".", alpha=0.5) + scale_fill_manual(values=c("4kHz_sup_v4"="#F89C74","5kHz_sup_v4"="#9EB9F3","5kHz_sup_v5"="#FE88B1"), labels=c("4kHz_sup_v4"="Dorado_4kHz_v4","5kHz_sup_v4"="Dorado_5kHz_v4","5kHz_sup_v5"="Dorado_5kHz_v5")) + scale_colour_manual(values=c("4kHz_sup_v4"="#F89C74","5kHz_sup_v4"="#9EB9F3","5kHz_sup_v5"="#FE88B1"), labels=c("4kHz_sup_v4"="Dorado_4kHz_v4","5kHz_sup_v4"="Dorado_5kHz_v4","5kHz_sup_v5"="Dorado_5kHz_v5")) + scale_y_continuous(limits=c(-50,100)) + xlab("Distance of C from CpG") + ylab("Difference in % methylation") + labs(title="", fill="Model", colour="Model") + my_theme + theme(axis.text.y=element_blank())

}


######################
# Combining the plots
######################

# Combining the plots profiling the Cs neighboring fully methylated CpGs together
plot_list[["us_100"]] + plot_list[["ds_100"]] + plot_layout(guides="collect", axis_titles="collect")

# Saving the plots as a png file
showtext_opts(dpi = 300)
ggsave('human_neighbouring_mod_jitter_box_CpG100.png', dpi = 300, bg = "white", height = 4, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()

# Combining the plots profiling the Cs neighboring fully unmethylated CpGs together
plot_list[["us_0"]] + plot_list[["ds_0"]] + plot_layout(guides="collect", axis_titles="collect")

# Saving the plots as a png file
showtext_opts(dpi = 300)
ggsave('human_neighbouring_mod_jitter_box_CpG0.png', dpi = 300, bg = "white", height = 4, width = 10)  
showtext_opts(dpi = 96)
showtext_auto()


