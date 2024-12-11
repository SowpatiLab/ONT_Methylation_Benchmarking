library(polars)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(showtext)
library(purrr)

setwd('') # set working directory
binDataSet <- function(fdf, xlabel, ylabel, bins=50){
    bins = as.numeric(bins)
    binsize <- 100/bins
    ef <- fdf$with_columns(
        c(
        (((pl$col(xlabel)/binsize)$floor()*binsize))$alias('x_bin'), # binning x axis to and flooring
        (((pl$col(ylabel)/binsize)$floor()*binsize))$alias('y_bin'), # binning y axis to and flooring
        pl$lit(1)$alias('count') # add counter column for aggregatuin sum
        )
    )$
        group_by(c('x_bin', 'y_bin', 'exp'))$   # group df by x and y labels
        agg(pl$col(c('count'))$sum())           # aggregation and sum
}
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color_palette <- rf(32)

# --------------------------- manual data curation

getTitles <- function(dataset, correl){
    select_exps <- dataset$exp %>% unique()
    correlframe <- as.data.frame(correl %>% filter( exp %in% select_exps))
    ret <- (
        correlframe %>% arrange(exp) 
        %>% mutate(titles=paste(exp, '\n r = ' , format(round(r, 3), nsmall = 3), ' \n n = ', format(n,  big.mark=",")))
        %>% select(c('exp', 'titles'))
    )
    namedVector <- setNames(ret$titles, ret$exp)
    return(namedVector)
}

plotBinnedHeatMap <- function(dataframeob, title=sprintf("%s contexts | bisulfite >= 1%%", SET_CONTEXT)){
    custom_labels <- getTitles(dataframeob, corr_df)
    
    p <- ggplot(dataframeob, aes(x=x_bin, y=y_bin, fill=count)) +
        geom_tile() + 
        scale_fill_gradientn(colors=color_palette, trans="log10") +
        facet_wrap(
        ~exp, ncol=length(unique(dataframeob$exp)),
        # labeller=c(1,2,3,4,5)
        labeller=labeller(exp=custom_labels)
        ) +
        labs(
        title=title,
        y="ONT Methylation", x="Bisulfite Methylation"
        )+
        theme(
        text = element_text(family = "lato", size = 20),
        plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 15),
        axis.text = element_text(size = 14.5),
        axis.title.x = element_text(vjust = 0),
        axis.title.y = element_text(hjust = 0.5),
        axis.title = element_text(size=24),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 20),
        panel.grid.major = element_line(colour="gray95"),
        strip.background = element_rect(fill="transparent")
        ) + coord_fixed(1)
    return(p)
}

# nonCpG context ===========================================================================================>

SET_CONTEXT='nonCpG' # context can be set to nonCpG, CHH and CHG
binned_full  <- as.data.frame(read.csv(sprintf('consolidated_binned/%s_binned.tsv', SET_CONTEXT), sep='\t')) # load connsolidated dataframe
corr_df      <- as.data.frame(read.csv(sprintf('consolidated_binned/%s_correl.tsv', SET_CONTEXT), sep='\t')) # 

head(binned_full)
head(corr_df)

showtext_opts(dpi = 300)

plotBinnedHeatMap(binned_full, title=sprintf("%s contexts | bisulfite full", SET_CONTEXT))
ggsave(sprintf('plots/faceted_correl_heatmap_%s_full.png',SET_CONTEXT), dpi = 300, bg = "white", height = 5, width =length(unique(binned_full$exp))*3)


