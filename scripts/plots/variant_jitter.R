# Plot Variant Calls
#
# Create a jitter plot for all variant calls in each pipeline.

library(ggplot2)


# Get arguments
args = commandArgs(trailingOnly=TRUE)

pipeline.name = args[1]
table.name = args[2]
out.jitter = args[3]
out.jitter.multi = args[4]
pipeline.name = args[5]


# Set pipeline name for plot titles
if (pipeline.name == "kestrel") {
    pipeline.name = "Kestrel"

} else if (pipeline.name == "gatk") {
    pipeline.name = "GATK"

} else {
    stop(paste("Unknown pipeline name:", pipeline.name))
}

# Read
df <- read.table(table.name, sep="\t", head=1, stringsAsFactors=FALSE)

# Make factors from call types (define order)
df$CALL = factor(df$CALL, levels=c("TP", "FP", "FN"))

# Get counts for greater than 30
df.hc <- subset(df, DEPTH > 30)
count.hc <- format(sapply(c("TP", "FP", "FN"), function(call.type) {sum(df.hc$CALL == call.type)}), big.mark=",", trim=TRUE)

# Get base plot
plot.base <- ggplot(df, aes(CALL, DEPTH)) +
    geom_jitter() +
    ylab("Alignment Depth at Variant Locus (reads)") +
    xlab(pipeline.name) +
    theme(panel.background=element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), axis.ticks.y=element_blank())
    scale_y_continuous(limits=c(-10, 600))


# Write stand-alone
plot.jitter <- plot.base +
    ggtitle(paste(pipeline.name, "Variant Calls")) +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(axis.text.x=element_text(size=16, color="gray30"), axis.text.y=element_text(size=12, color="black"))

ggsave(out.jitter, plot.jitter, height=5.5, width=5.5, units="in")

# Write multi-part
plot.jitter <- plot.base +
    theme(axis.text.x=element_text(size=18, color="black"), axis.text.y=element_text(size=14, color="black")) +
    theme(axis.title.x=element_text(size=18, color="gray30"), axis.title.y=element_text(size=18, color="black"))

ggsave(out.jitter.multi, plot.jitter, height=5.5, width=5.5, units="in")
