#############################
# Reading the input bed file 
#############################

data <- "exon_intersect.bed"	# This bed file contains methylation information for all seven tool-model variations at the corresponding exonic sites
df<-read.csv(data, sep="\t", header=F)
View(df)


#################################
# Loading the required libraries
#################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(showtext)
library(ggtext)
library(patchwork)

# Setting breaks for the colour gradient and setting the colour palette
my_breaks<-c(10^0,10^2,10^4)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral'))) # extend spectral color palette
r  <- rf(32)

# Rounding of the decimal numbers and counting the occurence for the heatmap binning
df_round<-df %>% mutate_if(is.numeric, round)
View(df_round)
names(df_round)<-c("chr","p1","p2",".","mod","strand","bis","sup_4kHz_v4","hac_5kHz_v4","sup_5kHz_v4","hac_5kHz_v5","sup_5kHz_v5","deepmod_4kHz","deepmod_5kHz")
df_melt<-melt(df_round, id.vars=c("chr","p1","p2",".","mod","strand","bis"))
View(df_melt)
df_count<-add_count(df_melt, bis, value)
View(df_count)


################################################################
# Calculating the correlation with the bisulfite reference data
################################################################

corr_4kHz_sup_v4<-cor(df$V7, df$V8, method="pearson")
print(corr_4kHz_sup_v4)
val1<-sprintf("%.3f", corr_4kHz_sup_v4)

corr_5kHz_sup_v4<-cor(df$V7, df$V10, method="pearson")
print(corr_5kHz_sup_v4)
val2<-sprintf("%.3f", corr_5kHz_sup_v4)

corr_5kHz_sup_v5<-cor(df$V7, df$V12, method="pearson")
print(corr_5kHz_sup_v5)
val3<-sprintf("%.3f", corr_5kHz_sup_v5)

corr_4kHz_deepmod<-cor(df$V7, df$V13, method="pearson")
print(corr_4kHz_deepmod)
val4<-sprintf("%.3f", corr_4kHz_deepmod)

corr_5kHz_deepmod<-cor(df$V7, df$V14, method="pearson")
print(corr_5kHz_deepmod)
val5<-sprintf("%.3f", corr_5kHz_deepmod)

corr_5kHz_hac_v4<-cor(df$V7, df$V9, method="pearson")
print(corr_5kHz_hac_v4)
val6<-sprintf("%.3f", corr_5kHz_hac_v4)

corr_5kHz_hac_v5<-cor(df$V7, df$V11, method="pearson")
print(corr_5kHz_hac_v5)
val7<-sprintf("%.3f", corr_5kHz_hac_v5)

# Setting the theme for the plots
my_theme <- theme(
	text = element_text(family = "lato", size = 14), 
	plot.title = element_text(size = 13.5, hjust=0.5), 
	axis.text = element_text(size = 12.5), 
	axis.title = element_text(size=13), 
	legend.text = element_text(size=12), 
	legend.title = element_text(size=14), 
	panel.grid.minor = element_blank(), 
	panel.background = element_blank(), 
	strip.text = element_text(size = 14),  
	panel.grid.major = element_line(colour="gray95")
)


#####################################################
# Generating the heatmap for each dataset separately
#####################################################

p1<-df_count %>% filter(variable=="deepmod_4kHz") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4), guide="none") + labs(title=paste("DeepMod2_4kHz_BiLSTM \n r=",val4, sep=""),y="ONT Methylation", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1) 

p2<-df_count %>% filter(variable=="sup_4kHz_v4") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4), guide="none") + labs(title=paste("Dorado_4kHz_v4sup \n r=",val1, sep=""),y="", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1) 

p3<-df_count %>% filter(variable=="deepmod_5kHz") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4)) + labs(title=paste("DeepMod2_5kHz_BiLSTM \n r=",val5, sep=""),y="", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1)  + theme(legend.position=c(2,0.5))

p4<-df_count %>% filter(variable=="sup_5kHz_v4") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4), guide="none") + labs(title=paste("Dorado_5kHz_v4sup \n r=",val2, sep=""),y="ONT Methylation", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1) 

p5<-df_count %>% filter(variable=="sup_5kHz_v5") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4), guide="none") + labs(title=paste("Dorado_5kHz_v5sup \n r=",val3, sep=""),y="", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1) 

p6<-df_count %>% filter(variable=="hac_5kHz_v4") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4), guide="none") + labs(title=paste("Dorado_5kHz_v4hac \n r=",val6, sep=""),y="", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1) 

p7<-df_count %>% filter(variable=="hac_5kHz_v5") %>% ggplot(aes(x=bis, y=value)) + scale_fill_gradientn(colors=r, trans="log10", breaks=my_breaks, limits=c(10^0,10^4), guide="none") + labs(title=paste("Dorado_5kHz_v5hac \n r=",val7, sep=""),y="", x="") + geom_bin_2d(binwidth = c(2,2)) + my_theme + coord_fixed(1) 


#######################################################################
# Setting the layout and combining all the heatmaps into a single plot
#######################################################################

layout=c(area(1,1),area(2,1),area(3,1),area(4,1))
plot_layout(design = layout)
patch1<-(p1 | p2 | p3 + plot_spacer()) + plot_layout()
patch4<-wrap_elements(patch1) +
    labs(tag = "Bisulfite Methylation") +
    theme(
        plot.tag.position = c(0.4,0.08)
    )

layout=c(area(2,1),area(2,2),area(2,3),area(2,4),area(2,5))
patch2<- (p4 | p6 | p5 | p7 + plot_spacer())
patch5<-wrap_elements(patch2) +
    labs(tag = "Bisulfite Methylation") +
    theme(
        plot.tag.position = c(0.4,0.08)
    )

layout=c(area(1,1,1,4),area(2,1,2,4))
plot_layout(design = layout)
patch3<- (patch4 / patch5) + plot_layout(guides="collect")   
patch3

# Saving the plot as a png file
showtext_opts(dpi = 300)
ggsave('exon_corr.png', dpi = 300, bg = "white", height = 6, width = 12)  
showtext_opts(dpi = 96)
showtext_auto()


