# Memory Usage Box Plot
#
# Create a box-plot of memory usage.

library(ggplot2)
library(reshape2)

# Get arguments
args = commandArgs(trailingOnly=TRUE)

filename.table <- args[1]
out.base.name <- args[2]

# Read data
df.size <- read.table(filename.table, header=TRUE, stringsAsFactors=FALSE)

df.size$seg_size <- NULL
df.size$seg_n <- NULL
df.size$bases <- NULL

df.size$reads = df.size$reads / 1e6
df.size$ikc_size = df.size$ikc_size / 1024
df.size$bam_size = df.size$bam_size / 1024

df.size <- melt(df.size, id.vars=c("accession", "reads"), measure_vars=c("ikc_size", "bam_size"))

# Get plot
plotObj <- ggplot(df.size, aes(reads, value, color=variable)) +
    xlab('Millions of 250 bp Reads') +
    ylab('Size of File (MB)') +
    geom_point(size=4, alpha=0.65) +
    scale_color_manual(values=c('bam_size'='turquoise3', 'ikc_size'='darkorchid3'), breaks=c('bam_size', 'ikc_size'), labels=c('BAM', 'IKC')) +
    guides(color = guide_legend(override.aes=list(alpha=1))) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    ggtitle('Size of BAM and IKC files by Millions of Reads') +
    theme(legend.justification=c(0, 1), legend.position=c(0, 1), legend.background=element_rect(color='black'), legend.key.size=unit(1, 'cm'), legend.title=element_blank()) +
    theme(axis.title=element_text(size=12, face='plain')) +
    theme(plot.title=element_text(size=16, hjust=0.5))

# Write
ggsave(paste(out.base.name, ".pdf", sep=""), plotObj, height=5.5, width=5.5, units="in")
ggsave(paste(out.base.name, ".eps", sep=""), plotObj, dev=cairo_ps, height=5.5, width=5.5, units="in")

save(plotObj, file=paste(out.base.name, ".RData", sep=""))
