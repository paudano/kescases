# Consensus size histogram
#
# Generate a histogram of the size of all consensus regions over all samples.

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

filename.con <- args[1]
filename.hap <- args[2]
filename.pdf <- args[3]

# Read table
df.con <- read.table(filename.con, header=TRUE)
df.hap <- read.table(filename.hap, header=TRUE)

df.con$Source <- "Consensus"
df.hap$Source <- "Haplotype"

df <- rbind(df.con, df.hap)

# Create plot
plotObj <- ggplot(df, aes(SIZE / 1e6, fill=Source)) +
    geom_histogram(bins=50, position="stack") +
    xlab("Consensus region size (Mbp)") +
    ylab("Number of samples") +
    ggtitle("Size of Consensus Regions") +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=10, color='black')) +
    theme(plot.title=element_text(size=24, hjust=0.5)) +
    scale_fill_manual(NULL, values=c('Consensus'='azure4', 'Haplotype'='seagreen3'))

# Write plot
ggsave(filename.pdf, plotObj, height=5, width=5, units="in")
