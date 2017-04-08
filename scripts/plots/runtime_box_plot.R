# Memory Usage Box Plot
#
# Create a box-plot of memory usage.

library(ggplot2)

# Get arguments
args = commandArgs(trailingOnly=TRUE)

table.kes.name = args[1]
table.gatk.name = args[2]
table.asm.name = args[3]
out.all = args[4]
out.noasm = args[4]

# Read tables
df.rt.kes <- read.table(table.kes.name, header=TRUE)
df.rt.gatk <- read.table(table.gatk.name, header=TRUE)
df.rt.asm <- read.table(table.asm.name, header=TRUE)

df.rt.kes$pipeline <- "Kestrel"
df.rt.gatk$pipeline <- "GATK"
df.rt.asm$pipeline <- "Assemble"

# Set table data
df.rt <- rbind(df.rt.kes, df.rt.gatk, df.rt.asm)
df.rt$pipeline <- factor(df.rt$pipeline, levels=c("Kestrel", "GATK", "Assemble"))
df.rt$user <- df.rt$user / 60
df.rt$real <- df.rt$real / 60

df.rt.noasm <- subset(df.rt, pipeline != "Assemble")

# Create plot - CPU - All
plot.rt <- ggplot(df.rt, aes(pipeline, user)) +
    geom_boxplot() +
    ggtitle("CPU Time by Pipeline") +
    ylab("CPU Time (m)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black"))

ggsave(out.all, plot.rt, height=5.5, width=5.5, units="in")

# Create plot - CPU - Kestrel vs GATK
plot.rt.noasm <- ggplot(df.rt.noasm, aes(pipeline, user)) +
    geom_boxplot() +
    ggtitle("CPU Time by Pipeline") +
    ylab("CPU Time (m)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black"))

ggsave(out.noasm, plot.rt.noasm, height=5.5, width=5.5, units="in")
