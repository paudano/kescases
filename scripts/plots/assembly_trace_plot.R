library(ggplot2)
library(tools)

# Get plot functions
source("scripts/plots/traceproc.R")

# Get arguments
args = commandArgs(trailingOnly=TRUE)

in.file.name = args[1]
out.file.name = args[2]

# Get plot
plotObj <- getTracePlot(in.file.name)

# Write
ggsave(out.file.name, plotObj, dev=cairo_ps, height=5.5, width=5.5, units="in")
