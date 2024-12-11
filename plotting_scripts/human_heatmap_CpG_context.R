library(polars)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(showtext)

my_breaks<-c(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

bin_size <- 50
context <- 'CpG'   # methylation context | can be set to CpG, CHH or CHG
setwd('') # set working directory to the appropriate location

loadDataSet <- function(sampling_rate, accuracy, context, model_version, xlabel, ylabel, bins=50){
    expNam <- sprintf('%skHz %s v%s', sampling_rate, accuracy, model_version)
    
    fdf <- pl$read_csv('data/CpG_5mC_20x_DoradoDeepmod2.bed', separator='\t')
    if (ylabel=="per"){
        
        if(accuracy!="deepmod2"){
        fdf <- pl$read_csv(sprintf('data/%s_HG002_%s_%skHz_v%s_%s.bed', context, accuracy, sampling_rate, model_version,'5mC'), separator='\t', has_header=FALSE)
        my_title <- sprintf("%skHz %s 5mC %s v%s", sampling_rate, accuracy, context, model_version)
        output_image = sprintf('correl_heatmap_%skHz_%s_v%s.png', sampling_rate, accuracy, model_version)
        } else {
        fdf <- pl$read_csv(sprintf('data/%s_HG002_%s_%skHz_%s.bed', context, accuracy, sampling_rate,'5mC'), separator='\t', has_header=FALSE)
        my_title <- sprintf("%skHz %s 5mC %s", sampling_rate, accuracy, context)
        output_image = sprintf('correl_heatmap_%skHz_%s.png', sampling_rate, accuracy)
        }
        
        fdf$columns <- c("chrom","p1","p2","mod","s","strand","cov","per","M","UM","per_bis","M_bis","UM_bis","context")
        fdf <- fdf$with_columns(c(
        (pl$col('M_bis') + pl$col('UM_bis'))$alias('cov_bis')
        ))$filter(pl$col('cov')>=20)$filter(pl$col('cov_bis')>=20)
        output_image = sprintf('correl_heatmap_%skHz_%s_v%s.png', sampling_rate, accuracy, model_version)
    } else {
        print('loading from combined dataset')
        fdf <- pl$read_csv('data/CpG_5mC_20x_DoradoDeepmod2.bed', separator='\t', has_header=TRUE)
        fdf$columns <- c("chrom","p1","p2","s","mod","strand","per_bis","4kHz_sup_v4","5kHz_hac_v4","5kHz_sup_v4","5kHz_hac_v5","5kHz_sup_v5","4kHz_deepmod2","5kHz_deepmod2","context")
        output_image = sprintf('correl_heatmap_intersected_%s.png', ylabel)
    }
    fdf <- fdf$select( pl$col( c(xlabel, ylabel)) )
    fdf <- fdf$with_columns( pl$lit(expNam)$alias('exp'))
    fdf$columns <- c('per_bis', 'per', 'exp')
    return(fdf)
}

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
        group_by(pl$col(c('x_bin', 'y_bin', 'exp')))$  # group df by x and y labels
        agg(pl$col(c('count'))$sum())           # aggregation and sum
}

getLookUp <- function(v){
    r <- format(round(corr[corr$exp==v, ]$r, 3), nsmall = 3)
    n <- format(corr[corr$exp==v, ]$n, big.mark=",")
    ret <- sprintf('%s\nr = %s\nn = %s', v,r,n)
    return(ret)
}

# =======================================================================
#  PLOTTING
# =======================================================================

d1 <- loadDataSet(4, "sup", context, 4, "per_bis", "per", 50)
d2 <- loadDataSet(5, "sup", context, 4, "per_bis", "per", 50)
d3 <- loadDataSet(5, "sup", context, 5, "per_bis", "per", 50)
d4 <- loadDataSet(4, "deepmod2", context, 5, "per_bis", "per", 50)
d5 <- loadDataSet(5, "deepmod2", context, 5, "per_bis", "per", 50)
d <- pl$concat(c(d1,d2,d3,d4,d5))

rename_exps <- data.frame(
    exp = c('5kHz deepmod2 v5', '4kHz deepmod2 v5', '4kHz sup v4', '5kHz sup v5', '5kHz sup v4'),
    new_value = c('5kHz_DeepMod2', '4kHz_DeepMod2', '4kHz_Dorado_v4', '5kHz_Dorado_v5', '5kHz_Dorado_v4')
)

d_binned <- binDataSet(d, "per_bis", "per")
d_binned <- as.data.frame(d_binned) %>% inner_join(rename_exps, by = "exp") %>% select(c("x_bin","y_bin","count", "new_value")) %>% rename(exp=new_value)

e1 <- as.data.frame(d)

corr <- e1 %>% group_by(exp) %>% summarize(r=cor(per, per_bis, method='pearson'), n=n(), .groups = 'drop')
corr <- corr %>% inner_join(rename_exps, by='exp') %>% select(c('r', 'n', 'new_value')) %>% rename(exp=new_value)

(
    ggplot() 
    + geom_tile(data=d_binned, aes(x=x_bin, y=y_bin, fill=count))
    + scale_fill_gradientn(colors=r, trans="log10") 
    + facet_wrap(
        ~exp, ncol=5,
        labeller=labeller(exp=getLookUp)
    )
    + coord_fixed(1)
    + labs(
        title='ONT 5mC performance on 5mC calling at CpG context',
        y="ONT Methylation", x="Bisulfite Methylation",
    )
    + theme(
        text = element_text(family = "lato", size = 20),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 12),
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
    )
)

showtext_opts(dpi = 300)
ggsave('faceted_correl_heatmap.png', dpi = 300, bg = "white", height = 5, width = 15)

head(d1)

# =======================================================================
#  GENERATE HEATMAPS FOR COMMONLY AVAILABLE SITES
# =======================================================================

d1 <- loadDataSet(4, "sup", context, 4, "per_bis", "4kHz_sup_v4", 50)
d2 <- loadDataSet(5, "sup", context, 4, "per_bis", "5kHz_sup_v4", 50)
d3 <- loadDataSet(5, "sup", context, 5, "per_bis", "5kHz_sup_v5", 50)
d4 <- loadDataSet(4, "deepmod2", context, 5, "per_bis", "4kHz_deepmod2", 50)
d5 <- loadDataSet(5, "deepmod2", context, 5, "per_bis", "5kHz_deepmod2", 50)
d <- pl$concat(c(d1,d2,d3,d4,d5))

rename_exps <- data.frame(
    exp = c('5kHz deepmod2 v5', '4kHz deepmod2 v5', '4kHz sup v4', '5kHz sup v5', '5kHz sup v4'),
    new_value = c('5kHz_DeepMod2', '4kHz_DeepMod2', '4kHz_Dorado_v4', '5kHz_Dorado_v5', '5kHz_Dorado_v4')
)

d_binned <- binDataSet(d, "per_bis", "per")
d_binned <- as.data.frame(d_binned) %>% inner_join(rename_exps, by = "exp") %>% select(c("x_bin","y_bin","count", "new_value")) %>% rename(exp=new_value)

e1 <- as.data.frame(d)
head(e1)

corr <- e1 %>% group_by(exp) %>% summarize(r=cor(per, per_bis, method='pearson'), n=n(), .groups = 'drop')
corr <- corr %>% inner_join(rename_exps, by='exp') %>% select(c('r', 'new_value')) %>% rename(exp=new_value)

getLookUpWithoutNumber <- function(v){
    r <- format(round( corr[corr$exp==v, ]$r, 3), nsmall = 3)
    print(r)
    ret <- sprintf('%s\nr = %s', v,r)
    return(ret)
}

(
    ggplot() 
    + geom_tile(data=d_binned, aes(x=x_bin, y=y_bin, fill=count))
    + scale_fill_gradientn(colors=r, trans="log10") 
    + facet_wrap(
        ~exp, ncol=5,
    )
    + coord_fixed(1)
    + labs(
        y="ONT Methylation", x="Bisulfite Methylation",
    )
    + theme(
        text = element_text(family = "lato", size = 20),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 12),
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
    )
)

showtext_opts(dpi = 300)
ggsave('faceted_correl_heatmap_intersected.png', dpi = 300, bg = "white", height = 4, width = 14)

