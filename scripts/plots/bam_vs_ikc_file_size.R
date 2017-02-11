# BAM vs IKC file size
#
# Create a scatter plot of BAM vs IKC file sizes for each sample.

library(ggplot2)


# Get arguments
args = commandArgs(trailingOnly=TRUE)

filename.table = args[1]
filename.blacklist = args[2]
out.base.name = args[3]


# Read sizes
df.size <- read.table(filename.table, header=TRUE, stringsAsFactors=FALSE)

df.size$ikc_size <- df.size$ikc_size / 1024
df.size$bam_size <- df.size$bam_size / 1024

max.size <- max(df.size$bam_size, df.size$ikc_size)

# Mark blacklisted samples
df.remregion <- read.table(filename.blacklist, sep='\t', header=TRUE, stringsAsFactors=FALSE)
df.remregion <- subset(df.remregion, start == 1 & end >= 2200000)  # Only use if the whole sample was blacklisted

df.size$removed <- ifelse(df.size$accession %in% df.remregion$accession, TRUE, FALSE)

# Plot - By sample
plot.size <- ggplot(df.size, aes(bam_size, ikc_size, color=removed)) +
    geom_point(size=2, alpha=0.75) +
    scale_x_continuous('BAM Size (MB)', limits=c(0, max.size)) +
    scale_y_continuous('IKC Size (MB)', limits=c(0, max.size * 0.25)) +
    scale_color_manual(NULL, labels=c('FALSE' = 'Analyzed', 'TRUE' = 'Removed'), values=c('FALSE' = 'black', 'TRUE' = 'firebrick')) +
    coord_fixed(ratio=1) +
    guides(color = guide_legend(override.aes=list(alpha=1))) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    ggtitle("Comparison of BAM and IKC file sizes by sample") +
    theme(legend.justification=c(1, 1), legend.position=c(1, 0.95), legend.background=element_rect(color='black'), legend.key.size=unit(0.5, 'cm'), legend.text=element_text(size=8)) +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=10, color='black')) +
    theme(plot.title=element_text(size=12, hjust=0.5))

# Write
ggsave(paste(out.base.name, ".pdf", sep=""), plot.size, height=5.5, width=5.5, units="in")
ggsave(paste(out.base.name, ".eps", sep=""), plot.size, dev=cairo_ps, height=5.5, width=5.5, units="in")

save(plot.size, file=paste(out.base.name, ".RData", sep=""))
