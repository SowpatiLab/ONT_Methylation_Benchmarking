#########################
# Reading the input file
#########################

TSS<-"df_indexed_TSS_5kb_flank_intersect_CpG_DoradoDeepmod2.bed"
df_TSS<-read.csv(TSS, sep="\t", header=TRUE)
# The header is: 'chr','p1','p2','feature','.','strand','ID','mod','per_bis','4kHz_sup_v4','5kHz_hac_v4','5kHz_sup_v4','5kHz_hac_v5','5kHz_sup_v5','4kHz_deepmod2','5kHz_deepmod2','context','index'
View(df_TSS)


#################################
# Loading the required libraries
#################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggalt)

df_melt = melt(df_TSS, id.vars = c('chr','p1','p2','feature','.','strand','ID','mod','context','index'))	# converting the data into long format
View(df_melt)
df_melt$bins_50 <- cut(df_melt$index, breaks=seq(-5001, 5001, 50))		# setting the bin size to 50bp
avg2 <- aggregate(value ~ bins_50 + variable, data = df_melt, FUN = mean)	# calculating the average methylation value of each bin
avg2_df<-as.data.frame(avg2)
names(avg2_df)<-c("bin","data","avg_value")
View(avg2_df)
write.csv(avg2_df, "avg2_df.csv")	# writing the dataframe into a csv file

# The csv file was modified on command line using awk for ease of usage
# awk -F ',' '{print $2,$3,$4,$5}' avg2_df.csv |sed 's/"//g'|sed 's/ -/,-/g'|sed 's/ /\t/g'|sed 's/]./\t/g'|sed 's/(//g'|sed 's/bin\t/bin_start\tbin_end\t/g' |sed 's/,/\t/g' > avg_50bp_bin.tsv

# It was then read back onto Rstudio

re_50_data<-"/path/to/avg_50bp_bin.tsv"

re_50_df<-read.csv(re_50_data, sep="\t", header=TRUE)
View(re_50_df)

# Setting the theme for plotting 
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
	panel.grid.major = element_line(colour="gray95"), 
	axis.text.x=element_text(angle=45,hjust=0.9)
)


##################################
# Plotting the data in a TSS plot
##################################

plot<-re_50_df %>% filter(data==c("per_bis","X4kHz_sup_v4","X5kHz_sup_v4", "X5kHz_sup_v5","X4kHz_deepmod2","X5kHz_deepmod2")) %>% ggplot(aes(x=bin_start, y=avg_value, color=data)) + labs(title="CpG methylation around TSS sites", x="Distance from TSS", y="% Methylation", color="Model") + geom_xspline(spline_shape = 1.5) + scale_colour_manual(breaks=c("per_bis","X4kHz_deepmod2","X4kHz_sup_v4","X5kHz_deepmod2","X5kHz_sup_v4","X5kHz_sup_v5"),labels=c("Bisulfite","DeepMod2_4kHz_BiLSTM","Dorado_4kHz_v4","DeepMod2_5kHz_BiLSTM","Dorado_5kHz_v4","Dorado_5kHz_v5"), values=c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#9EB9F3","#FE88B1")) + my_theme
plot

# Saving the plot in a png file
showtext_opts(dpi = 300)
ggsave('human_TSS_plot.png', dpi = 300, bg = "white", height = 4, width = 6)  
showtext_opts(dpi = 96)
showtext_auto()

