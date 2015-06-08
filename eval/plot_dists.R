#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2)) # for melt
suppressMessages(library(RColorBrewer))

plot.score.distribution <- function(d,output.file) {
    pdf(output.plot)
    md <- melt(d, id=c("index", "type"))
    source("simple_theme.R")
    p <- (ggplot(md, aes(x=value, fill=factor(type))) +
          ##geom_histogram(binwidth=2, position='identity') +          
          ##facet_grid(variable ~ ., scales="free", space="free") +
          geom_histogram(binwidth=2, position='dodge') +
          facet_wrap(~variable, ncol=1, scales="free_x") + 
          xlab("Score") +
          ylab("Count") +
          scale_fill_discrete(name="Dataset",
                              breaks=c("rand", "read", "shuff"),
                              labels=c("random", "original", "shuffled")) +
          simple_gdocs_theme()
          );
    plot(p)
    dev.off()
}

###
### if run from command line
if(!interactive()) {
    ##
    ## command line argument parsing
    args<-commandArgs(TRUE)
    input.csv <- args[1]
    output.plot <- args[2]
    plot.type <- args[3]
    if(length(args) != 3 && length(args) != 4) {
        print("USAGE: input.csv output.pdf [scores]");
    }
    max_rate <- if(length(args) == 4) as.numeric(args[4]) else 1.0
    ##
    ## plot different types
    ##    
    if(plot.type == "scores") {
        d <- read.csv(input.csv, header=TRUE)
        plot.score.distribution(d)
    } else {
        print("INVALID type specified");
    }
}
