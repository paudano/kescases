# Plot Variant Calls
#
# Create a jitter plot for all variant calls in each pipeline.

library(ggplot2)


# Get arguments
args = commandArgs(trailingOnly=TRUE)

pipeline.name = args[1]
table.name = args[2]
out.base.name = args[3]


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

# Full
plot.full <- ggplot(df, aes(CALL, DEPTH)) +
    geom_jitter() +
    ggtitle(paste(pipeline.name, "Variant Calls")) +
    ylab("Alignment Depth at Variant Locus (reads)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black")) +
    scale_y_continuous(limits=c(-10, 600))

# LC
plot.lc <- ggplot(subset(df, DEPTH <= 30), aes(CALL, DEPTH)) +
    geom_jitter() +
    ggtitle(paste(pipeline.name, "Variant Calls (Depth <= 30)")) +
    ylab("Alignment Depth at Variant Locus (reads)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    annotate("text", x=c(1, 2, 3), y=32, label=count.hc) +
    geom_segment(aes(x=1, y=33, xend=1, yend=35), arrow=arrow(length = unit(0.03, "npc"))) +
    geom_segment(aes(x=2, y=33, xend=2, yend=35), arrow=arrow(length = unit(0.03, "npc"))) +
    geom_segment(aes(x=3, y=33, xend=3, yend=35), arrow=arrow(length = unit(0.03, "npc"))) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black"))

# Write plots
ggsave(paste(out.base.name, "_all.pdf", sep=""), plot.full, height=5.5, width=5.5, units="in")
ggsave(paste(out.base.name, "_all.eps", sep=""), plot.full, dev=cairo_ps, height=5.5, width=5.5, units="in")

ggsave(paste(out.base.name, "_lc.pdf", sep=""), plot.lc, height=5.5, width=5.5, units="in")
ggsave(paste(out.base.name, "_lc.eps", sep=""), plot.lc, dev=cairo_ps, height=5.5, width=5.5, units="in")

save(plot.full, plot.lc, file=paste(out.base.name, ".RData", sep=""))
