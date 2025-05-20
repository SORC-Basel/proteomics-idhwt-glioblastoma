


dta_tmt_short             <- dta_tmt_log[,c(-1,-2,-3)]
dta_tmt_PCA_raw           <- data.frame(t(dta_tmt_short[]))
rownames(dta_tmt_PCA_raw) <- gsub("_exp\\d+$", "", rownames(dta_tmt_PCA_raw))



# Since we transpose the matrix to have the Protein names as columns
# Column names cannot start with numbers and cannot contain ;
# Thus we need to add an X before protein names with Numbers and replace ; with .
tmt_lr_of_interest_2 <- ifelse(grepl("^\\d", tmt_lr_of_interest), paste0("X", tmt_lr_of_interest), tmt_lr_of_interest)
tmt_lr_of_interest_2 <- gsub(";", ".", tmt_lr_of_interest_2)


dta_tmt_PCA <- dta_tmt_PCA_raw[, colnames(dta_tmt_PCA_raw) %in% tmt_lr_of_interest_2]
dta_tmt_PCA$groups <-   c(rep("Low Grade Glioma",length(cols_IDHLG)),
                                  rep("IDH HG", length(cols_IDHHG)),
                                  rep("Low Grade Glioma",  length(cols_IDHTERT)),
                                  rep("Low Grade Glioma",  length(cols_IDHCODEL)),
                                  rep("GB PN", length(cols_RTKPN)),
                                  rep("GB CL", length(cols_RTKCL)),
                                  rep("GB MES", length(cols_MES)) )


dta_tmt_PCA$groups <- factor(dta_tmt_PCA$groups, 
                                       levels = c("Low Grade Glioma", "IDH HG", "GB PN", "GB CL", "GB MES"))




# https://distill.pub/2016/misread-tsne/
# http://bioconductor.org/books/3.15/OSCA.basic/dimensionality-reduction.html
library("Rtsne")
tsne_full <- Rtsne(dta_tmt_PCA[,!(colnames(dta_tmt_PCA) %in% "groups")], max_iter = 10000,
              perplexity = floor((nrow(dta_tmt_PCA) - 1) / 3))


tsne_full <- Rtsne(dta_tmt_PCA[,!(colnames(dta_tmt_PCA) %in% "groups")], max_iter = 5000,
                   perplexity = 20, eta = 10)

plot(tsne_full$Y, col = dta_tmt_PCA$groups, pch=19)







library(umap)
plot_umap <- function(x, labels,
                      main="UMAP visualization",
                      # colors=c("#ff7f00", "#e377c2", "#17becf", "green", "red"),
                      # colors=c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba"),
                      # colors=c("#d7191c", "#fdae61", "#e377c2", "#abdda4", "#2b83ba"),
                      # colors=c("#EE7733", "#0077BB", "#EE3377", "#009988", "#BBBBBB"),
                      colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
                      # colors=c("#00008B", "#8B0000", "#008B00", "#8B008B", "#8B8B00"),
                      pad=0.1, cex=1.2, pch=19, add=FALSE, legend.suffix="",
                      cex.main=1, cex.legend=0.85) {
    
    layout <- x
    if (is(x, "umap")) {
        layout <- x$layout
    } 
    
    xylim <- range(layout)
    xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
    if (!add) {
        par(mar=c(0.2,0.7,1.2,0.7), ps=10)
        plot(xylim, xylim, type="n", axes=F, frame=F)
        rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
    }
    points(layout[,1], layout[,2], col=colors[as.integer(labels)],
           cex=cex, pch=pch)
    # mtext(side=3, main, cex=cex.main)
    mtext(at = -2, main, cex=cex.main)
    
    labels.u <- unique(labels)
    legend.pos <- "topright"
    legend.text <- as.character(labels.u)
    if (add) {
        legend.pos <- "bottomleft"
        legend.text <- paste(as.character(labels.u), legend.suffix)
    }
    
    legend(legend.pos, legend=legend.text, inset=0.03,
           col=colors[as.integer(labels.u)],
           bty="n", pch=unique(pch), cex=cex.legend)
}




plot_umap_ggplot <- function(x, labels,
                             main="UMAP visualization",
                             colors=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
                             shapes=c(16, 17, 18, 19, 20), # Default shapes
                             size=1.2) {
    
    # Prepare the data
    if (is(x, "umap")) {
        layout <- x$layout
    } else {
        layout <- x
    }
    
    # Create a data frame
    df <- data.frame(X = layout[, 1], Y = layout[, 2], Label = as.factor(labels))
    
    # Check if shapes length is less than number of unique labels and repeat if necessary
    unique_labels <- length(unique(labels))
    if(length(shapes) < unique_labels) {
        shapes <- rep(shapes, length.out = unique_labels)
    }
    
    # Create the plot
    p <- ggplot(df, aes(x = X, y = Y, Group = Label, color = Label, shape = Label, fill = Label)) +
        geom_point(size = size) +
        scale_color_manual(values = colors) +
        scale_shape_manual(values = shapes[1:unique_labels]) +
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_blank()) +
        ggtitle(main) +
        # theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_rect(colour = "gray", fill = NA, linewidth=0.1)
              ) +  
        theme(
            legend.text = element_text(size = 14),
            plot.title = element_text(size = 18)  # center title
        )
    
    # Return the plot
    return(p)
}


# Exclude groups 
umap_full <- umap(dta_tmt_PCA[,!(colnames(dta_tmt_PCA) %in% "groups")], preserve.seed = FALSE)

umap_DEP <- plot_umap_ggplot(umap_full$layout, dta_tmt_PCA$groups, size = 3,
                             shapes=c(16, 17, 15, 3, 7))

