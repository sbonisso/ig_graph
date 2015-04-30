
suppressMessages(library(ggthemes))

simple_gdocs_theme <- function() {
    (theme_gdocs() +
     theme(plot.background=element_rect(color="white", fill="white")) +
     theme(axis.title.y = element_text(angle=90)) +
     theme(axis.title.y = element_text(size=18)) +
     theme(axis.title.x = element_text(size=18)) +
     theme(axis.text = element_text(size=16)) 
     );
}
