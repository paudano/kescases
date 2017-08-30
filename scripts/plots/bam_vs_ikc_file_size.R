# BAM vs IKC file size
#
# Create a scatter plot of BAM vs IKC file sizes for each sample.

library(ggplot2)


# Get arguments
args = commandArgs(trailingOnly=TRUE)

filename.table = args[1]
filename.blacklist = args[2]
filename.out = args[3]
filename.out.eps = args[4]
filename.multi.out = args[5]


# Read sizes
df.size <- read.table(filename.table, header=TRUE, stringsAsFactors=FALSE)

df.size$ikc_size <- df.size$ikc_size / 1024
df.size$bam_size <- df.size$bam_size / 1024

max.size <- max(df.size$bam_size, df.size$ikc_size)
max.size.bam <- max(df.size$bam_size)
max.size.ikc <- max(df.size$ikc_size)

# Mark blacklisted samples
df.remregion <- read.table(filename.blacklist, sep='\t', header=TRUE, stringsAsFactors=FALSE)
df.remregion <- subset(df.remregion, start == 1 & end >= 2200000)  # Only use if the whole sample was blacklisted

df.size$removed <- ifelse(df.size$accession %in% df.remregion$accession, TRUE, FALSE)

# Plot - By sample
plot.base <- ggplot(df.size, aes(bam_size, ikc_size, color=removed)) +
    geom_point(size=4, alpha=0.75) +
    scale_x_continuous('BAM Size (MB)', limits=c(0, max.size.bam)) +
    scale_y_continuous('IKC Size (MB)', limits=c(0, max.size.ikc)) +
    scale_color_manual(NULL, labels=c('FALSE' = 'Analyzed', 'TRUE' = 'Removed'), values=c('FALSE' = 'black', 'TRUE' = 'firebrick')) +
    coord_fixed(ratio=1) +
    guides(color = guide_legend(override.aes=list(alpha=1))) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(legend.position="bottom") +
    theme(legend.background=element_rect(color='black'), legend.key.size=unit(1, 'cm'), legend.text=element_text(size=18))
#    theme(legend.justification=c(1, 0), legend.position=c(1, 0.05), legend.background=element_rect(color='black'), legend.key.size=unit(1, 'cm'), legend.text=element_text(size=18))


# Write stand-alone
plot.size <- plot.base +
    ggtitle("Comparison of BAM and IKC file sizes by sample") +
    theme(plot.title=element_text(size=24, hjust=0.5)) +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=10, color='black'))

ggsave(filename.out, plot.size, height=4.167, width=14, units="in")
ggsave(filename.out.eps, plot.size, dev=cairo_ps, height=4.167, width=14, units="in")


# Write multipart
plot.size <- plot.base +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=10, color='black')) +
    theme(axis.text.x=element_text(size=18, color="black"), axis.text.y=element_text(size=14, color="black")) +
    theme(axis.title.y=element_text(size=18, color="black"), axis.title.x=element_text(size=18, color="black"))

ggsave(filename.multi.out, plot.size, height=4.167, width=14, units="in")
