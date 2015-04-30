#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2)) # for melt
suppressMessages(library(RColorBrewer))

plot.heatmap <- function(d,output.plot) {
    pdf(output.plot)
    source("simple_theme.R")
    p <- (ggplot(d, aes(x=factor(k_d),y=factor(k_vj),fill=total))+
          geom_tile()+
          ##xlab("D") +
          xlab(expression("k"["D"])) +
          ##ylab("V/J")  +
          ylab(expression("k"["VJ"])) +
          ##scale_fill_gradientn(colours=brewer.pal(7,"Set2"),
          scale_fill_gradientn(colours=brewer.pal(10,"RdYlBu"),
                               limits=c(0,1.0)) +
          scale_x_discrete(expand = c(0,0),
                           breaks=seq(from=5,to=25,by=5),labels=seq(from=5,to=25,by=5)) +
          scale_y_discrete(expand = c(0,0),
                           breaks=seq(from=5,to=25,by=5),labels=seq(from=5,to=25,by=5)) +
          simple_gdocs_theme());
    plot(p);
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
    if(length(args) != 3) {
        print("USAGE: input.csv output.pdf [heatmap]");
        stopifnot(length(args) == 4);
    }
    ##
    ## read data
    ##
    d <- read.csv(input.csv, sep=",", header=TRUE)  
    ##dM.dat <- d
    print(dim(d))
    print(names(d))
    ##
    ## plot the clone expression
    if(plot.type == "heatmap") {
        plot.heatmap(d,output.plot);
    } else {
        print("INVALID plot type specified!");
    }
}
