# Memory Usage Box Plot
#
# Create a box-plot of memory usage.

library(ggplot2)


# Get arguments
args = commandArgs(trailingOnly=TRUE)

table.kes.name = args[1]
table.gatk.name = args[2]
table.asm.name = args[3]
out.mem = args[4]
out.mem.multi = args[5]

# Read tables
df.mem.kes <- read.table(table.kes.name, header=TRUE)
df.mem.gatk <- read.table(table.gatk.name, header=TRUE)
df.mem.asm <- read.table(table.asm.name, header=TRUE)

df.mem.kes$pipeline <- "Kestrel"
df.mem.gatk$pipeline <- "GATK"
df.mem.asm$pipeline <- "Assemble"

df.mem <- rbind(df.mem.kes, df.mem.gatk, df.mem.asm)

df.mem$pipeline <- factor(df.mem$pipeline, levels=c("Kestrel", "GATK", "Assemble"))

# Transform axes
df.mem$rss = df.mem$rss / 1024^2 # KB to GB

# Create plot
plot.mem.base <- ggplot(df.mem, aes(pipeline, rss)) +
    geom_boxplot() +
    ylab("Maximum Resident Memory (GB)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black"))

# Save stand-alone
plot.mem <- plot.mem.base +
    ggtitle("Maximum Memory Usage by Pipeline") +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black")) +
    theme(plot.title=element_text(size=16, hjust=0.5))

ggsave(out.mem, plot.mem, height=5.5, width=5.5, units="in")

# Save multi-part
plot.mem <- plot.mem.base +
    theme(axis.text.x=element_text(size=18, color="black"), axis.text.y=element_text(size=14, color="black")) +
    theme(axis.title.y=element_text(size=18, color="black"))

ggsave(out.mem.multi, plot.mem, height=5.5, width=5.5, units="in")
