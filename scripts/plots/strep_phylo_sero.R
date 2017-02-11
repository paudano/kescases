library(circlize)

# Get arguments
args = commandArgs(trailingOnly=TRUE)

filename.ani <- args[1]
filename.sra <- args[2]
filename.out <- args[3]


# Get data for dendogram and heatmap
df.ani <- read.table(filename.ani, check.names=FALSE)

hc <- hclust(dist(1 - df.ani))
labels <- hc$labels
n <- length(labels)

dend <- as.dendrogram(hc)
max_height <- attr(dend, "height")

labels <- labels[order.dendrogram(dend)]

# Get colors for serotypes
df.sra <- read.table(filename.sra, sep='\t', header=TRUE, stringsAsFactors=FALSE)
df.sra <- df.sra[ ,c('Run_s', 'serotype_s')]
colnames(df.sra) <- c('accession', 'serotype')
rownames(df.sra) <- df.sra$accession

df.sra$serotypemajor <- sub('(\\d+).*', '\\1', df.sra$serotype)  # Remove sub-types

seronames <- sort(unique(df.sra$serotypemajor))

serocol <- rainbow(length(seronames), s=0.7, v=1.0)
names(serocol) <- seronames

# Add colors to df.sra
df.sra$col <- sapply(df.sra$serotypemajor, function(majtype){return(serocol[majtype])})

# Add reference
df.sra <- rbind(df.sra, data.frame(row.names='NC_003028', accession='NC_003028', serotype='REF', serotypemajor='NC_003028', col='#FFFFFF'))

# Init circos
pdf(filename.out, width=10, height=10)

par(mar = c(0, 0, 0, 0))
circos.clear()
circos.par("gap.degree" = 0, cell.padding=c(0.02, 0, 0.02, 0))
circos.initialize(factors='a', xlim=c(0, n))

# Serotype
circos.trackPlotRegion(ylim=c(0, 1), bg.border=NA, track.height=0.3, track.margin=c(0, 0.01), cell.padding=c(0, 0, 0.02, 0),
    panel.fun=function(x, y) {
        for (i in seq_len(n)) {
            circos.rect(i - 1, 0, i, 1, col=df.sra[labels[i],]$col)
            circos.text(i - 0.5, 0.8, df.sra[labels[i],]$serotype, cex=0.7, facing='clockwise', niceFacing=TRUE)
        }
    }
)

# Draw dendrogram
circos.trackPlotRegion(ylim=c(0, max_height), bg.border=NA, track.height=0.65, track.margin=c(0.01, 0), cell.padding=c(0.02, 0, 0, 0),
    panel.fun=function(x, y) {
        circos.dendrogram(dend, max_height=max_height)
    }
)

# Clear
circos.clear()
