library(ggplot2)

# getTracePlot
#
# Get the ggplot object of a memory trace (line plot, time on x-axis, memory usage on y-axis).
#
# Params
#   * file.name: Output from traceproc
#   * xlabel: X-axis label.
#   * ylabel: Y-axis label.
getTracePlot <- function(file.name, xlabel="Time (s)", ylabel="Memory (GB)", rotate.x=FALSE) {

    # Read
    df <- read.table(file.name, comment.char='', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

    # RSS kB to GB
    df["rss_kb"] <- df["rss_kb"] / 1024^2

    # Plot
    plotObj <- ggplot(df, aes(x=time_ms, y=rss_kb)) +
        geom_line() +
        xlab(xlabel) +
        ylab(ylabel) +
        scale_x_continuous(labels = scales::comma)

    # Rotate label
    if (rotate.x) {
        plotObj <- plotObj + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    return(plotObj)
}
