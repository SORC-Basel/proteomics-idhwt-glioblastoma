



dta_up_down_lr_venn <- cbind(dta_up_down_lr, !is.na(dta_up_down_lr[,2:5]))

list_lr <- list(
    "IDH HG vs Low grade glioma" = dta_up_down_lr_venn$tmt_lr_of_interest[dta_up_down_lr_venn[6] == TRUE],
    "GB PN vs Low grade glioma"  = dta_up_down_lr_venn$tmt_lr_of_interest[dta_up_down_lr_venn[7] == TRUE],
    "GB CL vs Low grade glioma"  = dta_up_down_lr_venn$tmt_lr_of_interest[dta_up_down_lr_venn[8] == TRUE],
    "GB MES vs Low grade glioma" = dta_up_down_lr_venn$tmt_lr_of_interest[dta_up_down_lr_venn[9] == TRUE]
)


# VennDiagram::venn.diagram(x = list_lr)
# gplots::venn(dta_up_down_lr_venn[,6:9])
p_venn <- ggVennDiagram::ggVennDiagram(list_lr, label_alpha = 0, label = "count", label_geom = "label", label_color = "black", edge_size = 0.5, 
                                       set_size   = 5)  + 
    scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
    scale_color_manual(values=c("#000000", "#000000", "#000000", "#000000")) +
    theme(legend.position = "none",
          plot.margin = margin(0, 10, 0, 10),
          # text = element_text(size = 90),
          plot.title = element_text(size = 20)) +
    labs(title = "Venn diagram - High grade glioma subgroups vs Low grade glioma") +  scale_x_continuous(expand = expansion(mult = 0.2)) + 
    theme(plot.title = element_text(hjust = 0.5))
# scale_fill_distiller(palette = "RdBu")



## Venn Diagramm only UP regulted
dta_up_lr_venn <- dta_up_down_lr_venn[apply(dta_up_down_lr_venn[,2:5] == "Up", 1, any),]
dta_up_lr_venn <- dta_up_lr_venn[!is.na(dta_up_lr_venn$tmt_lr_of_interest ), ]

list_lr_up <- list(
    "IDH HG vs Low grade glioma" = dta_up_lr_venn$tmt_lr_of_interest[dta_up_lr_venn[6] == TRUE],
    "GB PN vs Low grade glioma"  = dta_up_lr_venn$tmt_lr_of_interest[dta_up_lr_venn[7] == TRUE],
    "GB CL vs Low grade glioma"  = dta_up_lr_venn$tmt_lr_of_interest[dta_up_lr_venn[8] == TRUE],
    "GB MES vs Low grade glioma" = dta_up_lr_venn$tmt_lr_of_interest[dta_up_lr_venn[9] == TRUE]
)


p_venn_up <- ggVennDiagram::ggVennDiagram(list_lr_up, label_alpha = 0, label = "count", label_geom = "label", label_color = "black", edge_size = 0.5)  + 
    scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
    scale_color_manual(values=c("#000000", "#000000", "#000000", "#000000")) +
    theme(legend.position = "none") +
    labs(title = "Venn diagram - High grade glioma subgroups vs Low grade glioma",
         subtitle = "Only up-regulated proteins") +  scale_x_continuous(expand = expansion(mult = 0.2)) + 
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5)
          )
# scale_fill_distiller(palette = "RdBu")
# scale_fill_distiller(palette = "RdBu")




## Venn Diagramm only DOWN regulted
dta_down_lr_venn <- dta_up_down_lr_venn[apply(dta_up_down_lr_venn[,2:5] == "Down", 1, any),]
dta_down_lr_venn <- dta_down_lr_venn[!is.na(dta_down_lr_venn$tmt_lr_of_interest ), ]

list_lr_down <- list(
    "IDH HG vs Low grade glioma" = dta_down_lr_venn$tmt_lr_of_interest[dta_down_lr_venn[6] == TRUE],
    "GB PN vs Low grade glioma"  = dta_down_lr_venn$tmt_lr_of_interest[dta_down_lr_venn[7] == TRUE],
    "GB CL vs Low grade glioma"  = dta_down_lr_venn$tmt_lr_of_interest[dta_down_lr_venn[8] == TRUE],
    "GB MES vs Low grade glioma" = dta_down_lr_venn$tmt_lr_of_interest[dta_down_lr_venn[9] == TRUE]
)


p_venn_down <- ggVennDiagram::ggVennDiagram(list_lr_down, label_alpha = 0, label = "count", label_geom = "label", label_color = "black", edge_size = 0.5)  + 
    scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") +
    scale_color_manual(values=c("#000000", "#000000", "#000000", "#000000")) +
    theme(legend.position = "none") +
    labs(title = "Venn diagram - High grade glioma subgroups vs Low grade glioma",
         subtitle = "Only down-regulated proteins") +  scale_x_continuous(expand = expansion(mult = 0.2)) + 
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5)
    )
# scale_fill_distiller(palette = "RdBu")








up_or_down <- as.data.frame(dta_up_down_lr_venn[,1:5])


# Define a function to evaluate each row and determine if it's Up, Down, or Mixed
evaluate_row <- function(row) {
    # Get the non-NA values from the row
    non_na_values <- row[!is.na(row)]
    
    # Count the number of "Up" and "Down" values
    up_count <- sum(non_na_values == "Up")
    down_count <- sum(non_na_values == "Down")
    
    # If there are no non-NA values, return NA
    if(length(non_na_values) == 0) return(c("NA", 0, 0))
    
    # Check if all non-NA values are "Up"
    if(all(non_na_values == "Up")) return(c("Up", up_count, down_count))
    
    # Check if all non-NA values are "Down"
    if(all(non_na_values == "Down")) return(c("Down", up_count, down_count))
    
    # If there is a mix of "Up" and "Down", return "Mixed"
    return(c("Mixed", up_count, down_count))
}
# Apply the function to each row and create a new column "Direction" with the results
# up_or_down$Direction <- apply(up_or_down[-1], 1, evaluate_row)
# table(up_or_down$Direction, useNA = "always")

results <- t(as.data.frame(apply(up_or_down[-1], 1, evaluate_row)))
colnames(results) <- c("Direction", "Up_Count", "Down_Count")

# Bind the new columns to the original data frame
up_or_down <- cbind(up_or_down, results)



