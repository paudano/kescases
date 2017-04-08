# Memory Usage Box Plot
#
# Create a box-plot of memory usage.

library(ggplot2)


# Get arguments
args = commandArgs(trailingOnly=TRUE)

table.kes.name = args[1]
table.gatk.name = args[2]
table.asm.name = args[3]
out.name = args[4]

# Read tables
df.mem.kes <- read.table(table.kes.name, header=TRUE)
df.mem.gatk <- read.table(table.gatk.name, header=TRUE)
df.mem.asm <- read.table(table.asm.name, header=TRUE)

df.mem.kes$pipeline <- "Kestrel"
df.mem.gatk$pipeline <- "GATK"
df.mem.asm$pipeline <- "Assemble"

df.mem <- rbind(df.mem.kes, df.mem.gatk, df.mem.asm)

# Transform axes
df.mem$rss = df.mem$rss / 1024^2 # KB to GB

# Create plot
plot.mem <- ggplot(df.mem, aes(pipeline, rss)) +
    geom_boxplot() +
    ggtitle("Maximum Memory Usage by Pipeline") +
    ylab("Maximum Resident Memory (GB)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black"))

# Write
ggsave(out.name, plot.mem, height=5.5, width=5.5, units="in")
