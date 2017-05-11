# Plot variant depth (by caller) vs alignment depth at the variant locus.

library(ggplot2)


# Get arguments
args <- commandArgs(trailingOnly=TRUE)

table.name <- args[1]
out.pdf <- args[2]
pipeline <- args[3]

if (pipeline == 'kestrel') {
    pipeline <- 'Kestrel'
} else if (pipeline == 'gatk') {
    pipeline <- 'GATK'
}

# Read tables
df.var <- read.table(table.name, header=TRUE, sep="\t")

# Create plot
plot.depth <- ggplot(df.var, aes(x = DP, y = DEPTH)) +
    geom_point() +
    ggtitle(paste("Call vs Alignment Depth -", pipeline)) +
    xlab("Variant Call Depth") +
    ylab("Alignment Depth") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=10, color='black'))

# Write
ggsave(out.pdf, plot.depth, height=5.5, width=5.5, units="in")
