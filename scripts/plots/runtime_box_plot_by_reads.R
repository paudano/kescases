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
out.noasm = args[5]
out.all.multi = args[6]
out.noasm.multi = args[7]

# Read tables
df.rt.kes <- read.table(table.kes.name, header=TRUE)
df.rt.gatk <- read.table(table.gatk.name, header=TRUE)
df.rt.asm <- read.table(table.asm.name, header=TRUE)

df.rt.kes$pipeline <- "Kestrel"
df.rt.gatk$pipeline <- "GATK"
df.rt.asm$pipeline <- "Assemble"

# Divide by millions of reads
df.rt.kes$user = df.rt.kes$user / (df.rt.kes$reads / 1e6)
df.rt.gatk$user = df.rt.gatk$user / (df.rt.gatk$reads / 1e6)
df.rt.asm$user = df.rt.asm$user / (df.rt.asm$reads / 1e6)

# Set table data
df.rt <- rbind(df.rt.kes, df.rt.gatk, df.rt.asm)
df.rt$pipeline <- factor(df.rt$pipeline, levels=c("Kestrel", "GATK", "Assemble"))
df.rt$user <- df.rt$user / 60
df.rt$real <- df.rt$real / 60

df.rt.noasm <- subset(df.rt, pipeline != "Assemble")


### Create plot - CPU - All ###
plot.rt.base <- ggplot(df.rt, aes(pipeline, user)) +
    geom_boxplot() +
    ylab("CPU Time per Million Reads (m)") +
    xlab("") +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())


# Save stand-alone
plot.rt <- plot.rt.base +
    ggtitle("CPU Time by Pipeline") +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black")) +
    theme(plot.title=element_text(size=16, hjust=0.5))

ggsave(out.all, plot.rt, height=5.5, width=5.5, units="in")

# Save multi-part
plot.rt <- plot.rt.base +
    theme(axis.text.x=element_text(size=18, color="black"), axis.text.y=element_text(size=14, color="black")) +
    theme(axis.title.y=element_text(size=18, color="black"))

ggsave(out.all.multi, plot.rt, height=5.5, width=5.5, units="in")


### Create plot - CPU - Kestrel vs GATK ###
plot.rt.base <- ggplot(df.rt.noasm, aes(pipeline, user)) +
    geom_boxplot() +
    ggtitle("CPU Time by Pipeline") +
    ylab("CPU Time per Million Reads (m)") +
    xlab("") +
    theme(plot.title=element_text(size=16, hjust=0.5)) +
    theme(panel.background = element_blank(), axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(panel.grid.major.y=element_line(color="gray"), panel.grid.minor.y=element_line(color="gray"), axis.ticks.y=element_blank()) +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black"))

# Save stand-alone
plot.rt <- plot.rt.base +
    ggtitle("CPU Time by Pipeline") +
    theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=12, color="black")) +
    theme(plot.title=element_text(size=16, hjust=0.5))

ggsave(out.noasm, plot.rt, height=5.5, width=5.5, units="in")

# Save multi-part
plot.rt <- plot.rt.base +
    theme(axis.text.x=element_text(size=18, color="black"), axis.text.y=element_text(size=14, color="black")) +
    theme(axis.title.y=element_text(size=18, color="black"))

ggsave(out.noasm.multi, plot.rt, height=5.5, width=5.5, units="in")
