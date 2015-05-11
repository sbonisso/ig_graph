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
plot.full.heatmap <- function(d,output.plot) {
    pdf(output.plot)
    source("simple_theme.R")
    d$X.name <- factor(d$X.name, levels(d$X.name)[order(d$X.name)])
    d$variable <- factor(d$variable, levels(d$variable)[order(d$X.name)])
    p <- (ggplot(d, aes(x=X.name,y=variable,fill=value))+
          geom_tile()+
          geom_text(aes(X.name, variable, label = value), color = "#073642", size = 4) +
          ylab("") + xlab("") + 
          ## scale_fill_gradientn(colours=brewer.pal(10,"RdYlBu"),
          ##                       limits=c(0,1.0)) +
          scale_fill_gradient(name="Jaccard index",
                              low = "#fdf6e3", high = "steelblue",
                              breaks=seq(0, 1, by = 0.2), limits = c(0.3, 1)) +
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +          
          ##simple_gdocs_theme() +
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top",
                     title.hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.ticks = element_blank(),
                legend.justification = c(1, 0),
                ##legend.position = c(0.9, 1.05),
                legend.position = "top",
                legend.direction = "horizontal") +
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", 
                     title.hjust = 0.5))
          );
    plot(p);
    dev.off()
}
plot.half.heatmap <- function(d,output.plot) {
    pdf(output.plot)
    source("simple_theme.R")
    p <- (ggplot(d, aes(x=X.name,y=variable,fill=value))+
          geom_tile()+
          geom_text(aes(X.name, variable, label = value), color = "#073642", size = 4) +
          ylab("") + xlab("") + 
          scale_fill_gradient(name="Jaccard index",
                              low = "#fdf6e3", high = "steelblue",
                              breaks=seq(0, 1, by = 0.2), limits = c(0.2, 1)) +
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +          
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top",
                     title.hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.ticks = element_blank(),
                legend.justification = c(1, 0),
                          legend.position = "top",
                legend.direction = "horizontal") +
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", 
                     title.hjust = 0.5))
          );
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
    if(length(args) != 3 && length(args) != 4) {
        print("USAGE: input.csv output.pdf [heatmap/half_mat] [max_rate]");
        ##stopifnot(length(args) == 4);
    }
    max_rate <- if(length(args) == 4) as.numeric(args[4]) else 1.0
    ##
    ## read data
    ##
    
    ##
    ## plot the clone expression
    if(plot.type == "heatmap") {
        d <- read.csv(input.csv, sep=",", header=TRUE)  
        d$total <- d$total / max_rate;
        plot.heatmap(d,output.plot);
    } else if(plot.type == "half_mat") {
        ## read in as matrix
        d <- read.csv(input.csv, sep=",", header=TRUE, row.names=1)
        d1 <- as.matrix(d)
        ##print(d1);
        ## remove lower triangle, add in diagonal
        d1[lower.tri(d1)] <- NA
        diag(d1) <- 1;
        d2 <- data.frame(d1)
        d2$X.name <- rownames(d2);
        ##print(d2);
        ## melt to long format
        md <- melt(d2);
        ##print(md)
        ## remove NA 
        md <- md[-which(is.na(md$value)),]
        md2 <- data.frame(md)
        
        ## sort names by frequency
        tn <- sort(table(md2$X.name), decreasing=TRUE)
        
        ## set as factor levels - based on freq
        md2$X.name <- factor(md2$X.name, levels=factor(names(tn)))
        md2$variable <- factor(md2$variable, levels=factor(names(tn)))
        ##md2$X.name <- factor(md2$X.name, levels=c("VDJsolver","AbOrigin","IgBLAST","IMGT", "iHMMune", "JOINSOLVER"))
        ##md2$variable <- factor(md2$variable, levels=c("VDJsolver","AbOrigin","IgBLAST","IMGT", "iHMMune", "JOINSOLVER"))
        ##print(levels(md2$X.name))
        ##print(levels(md2$variable))
        
        plot.half.heatmap(md2,output.plot);
    } else {
        print("INVALID plot type specified!");
    }
}
