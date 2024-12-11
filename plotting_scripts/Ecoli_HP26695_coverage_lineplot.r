########################
# Setting the variables
########################

mod <- '6mA'	#set to 5mC or 6mA
motif <- 'ATTAAT'	#set to the motif being profiled
species <- 'H.pylori'	# set to E.coli or H.pylori

#########################
# Reading the input file
#########################

data<-sprintf('metrics_%s_%s_%s.csv', motif, mod, species)

df<-read.csv(data, header=T)
View(df)
# The header has the following columns
# True negative, False positive, True positive, False negative, Precision, Recall, F1, Specificity, Set, Coverage, Data, Mode, Model, Context, Motif, AP, AN, Accuracy


#################################
# Loading the required libraries
#################################

library(ggplot2)
library(ggtext)
library(showtext)

######################################################################################
# Function to calculate the standard deviation around the mean for error bar plotting
######################################################################################

data_summary <- function(data, varname, groupnames){
 	require(plyr)
 	summary_func <- function(x, col){
 	    c(mean = mean(x[[col]], na.rm=TRUE),
 	      sd = sd(x[[col]], na.rm=TRUE))
 	}
 	data_sum<-ddply(data, groupnames, .fun=summary_func,
 	                varname)
 	data_sum <- rename(data_sum, c("mean" = varname))
 	return(data_sum)
}

# Setting the theme, title and the output file names
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

my_title <- sprintf("Species: %s | Mod:%s | Context:%s", species, mod, motif)

##########################################################################################
# Looping through the performance metrics and plotting the respective coverage line plots
##########################################################################################

array <- c("F1", "Precision", "Recall", "Specificity")

for (metric in array) {

	# Generating the dataframes with the summary function information for each of the performance metrics
	df_metric <- data_summary(df, varname=paste(metric), groupnames=c("Data", "Mode","Model","Coverage"))
	View(df_metric)


	out_img <- sprintf("%s_%s_%s_%s_coverage_lineplot.png", species, mod, motif, metric)


	# Plotting the data
	plot<-ggplot(df_metric, aes(x=Coverage, y=.data[[metric]], group=interaction(Data,Mode,Model), color=interaction(Data,Mode,Model))) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=.data[[metric]]-sd, ymax=.data[[metric]]+sd), width=.2, position=position_dodge(0.05)) + scale_x_discrete(limits=c("5x","10x","20x","30x","40x","50x","60x","70x")) + scale_colour_manual(values=c("#F89C74","#9EB9F3","#FE88B1"), labels=c("4kHzsup_v4","5kHzsup_v4","5kHzsup_v5")) + labs(title=my_title, x="Coverage", y=paste(metric), color="Model") + my_theme
	plot
	
	# Saving the plot as a png file
	showtext_opts(dpi = 300)
	ggsave(out_img, dpi = 300, bg = "white", height = 3, width = 6)  
	showtext_opts(dpi = 96)	
	showtext_auto()
}

