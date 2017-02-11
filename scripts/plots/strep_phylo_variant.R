library(circlize)

# Get arguments
args = commandArgs(trailingOnly=TRUE)

filename.ani <- args[1]
filename.var <- args[2]
filename.out <- args[3]

# Get data for dendogram and heatmap
df.ani <- read.table(filename.ani, check.names=FALSE)

hc <- hclust(dist(1 - df.ani))
labels <- hc$labels
n <- length(labels)

dend <- as.dendrogram(hc)
max_height <- attr(dend, "height")

labels <- labels[order.dendrogram(dend)]

# Get data for variants
df.var <- read.table(filename.var, sep='\t', header=TRUE, stringsAsFactors=FALSE)
df.var <- subset(df.var, CALL == 'TP' | CALL == 'FN')

varcount <- table(df.var$SAMPLE)

# Set variant count to 0 for samples with no variants
for (i in seq_len(length(labels))) {

    if (! labels[i] %in% names(varcount)) {
        varcount[labels[i]] <- 0
    }
}

# Heatmap colors
colfunc <- colorRampPalette(c('red', 'yellow', 'blue'))
col <- colfunc(100)

sample.col <- as.numeric(df.ani['NC_003028',])
sample.col[sample.col == 1] <- NA
sample.col.min <- min(sample.col, na.rm=TRUE)
sample.col.max <- max(sample.col, na.rm=TRUE)
sample.col <- floor((sample.col  - sample.col.min) / (sample.col.max - sample.col.min) * 99) + 1
sample.col <- sapply(sample.col, function(val){col[val]})
names(sample.col) <- colnames(df.ani['NC_003028',])
sample.col['NC_003028'] <- '#FFFFFF'

# Init circos

pdf(filename.out, width=10, height=10)

par(mar = c(0, 0, 0, 0))
circos.clear()
circos.par("gap.degree" = 0, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors='a', xlim=c(0, n))

# Bar plot
circos.trackPlotRegion(ylim=c(0, max(varcount)), bg.border=NA, track.height=0.3,
    panel.fun=function(x, y) {
        for (i in seq_len(n)) {
            if (varcount[labels[i]] > 0) {
                circos.lines(i - 0.5, varcount[labels[i]], type='h', straight=TRUE, lwd=5)
            } else {
                circos.lines(i - 0.5, 0, type='h', straight=TRUE, lwd=5, col='gray')
            }
        }
    }
)

# Heatmap
circos.trackPlotRegion(ylim=c(0, 1), bg.border=NA, track.height=0.15, track.margin=c(0.00, 0.01),
    panel.fun=function(x, y) {
        for (i in seq_len(n)) {
            circos.rect(i - 1, 0, i, 1, col=sample.col[labels[i]])
        }
    }
)

# Draw dendrogram
circos.trackPlotRegion(ylim=c(0, max_height), bg.border=NA, track.height=0.45, track.margin=c(0.01, 0.00),
    panel.fun=function(x, y) {
        circos.dendrogram(dend, max_height=max_height)
    }
)

# Clear
circos.clear()

# Finish plot
dev.off()
