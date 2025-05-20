# ----------------------------------
# Heatmap - significant proteins
# ----------------------------------
## Heatmap of glioma subgroups 
## Colour palette
my_palette <- gplots::colorpanel(256, "blue2", "white", "red2")

tmt_group_means <- fit_tmt_lr_1_full$coefficients
tmt_group_means <- tmt_group_means[rownames(tmt_group_means) %in% tmt_lr_of_interest,]

# Manual scaling of values for heatmap
tmt_lr_means_scaled <- tmt_group_means[,c("LowGrade", "IDHHG", "RTKPN", "RTKCL",  "MES")]
tmt_lr_means_scaled = t(apply(tmt_lr_means_scaled, 1, scale))

# Order of Rows if dendofram is switched of
colnames(tmt_lr_means_scaled) <- c("Low Grade Glioma", "IDH HG", "GB PN", "GB CL", "GB MES")
rownames(tmt_lr_means_scaled) <- str_remove(rownames(tmt_lr_means_scaled), "_HUMAN")
rownames(tmt_lr_means_scaled) <- str_replace(rownames(tmt_lr_means_scaled), "VIME_CON-HUMAN,VIME", "VIME")
rownames(tmt_lr_means_scaled) <- str_replace(rownames(tmt_lr_means_scaled), "LV746;LV743_HUMAN", "IGLV746")


heatmap_DEP <- Heatmap(tmt_lr_means_scaled, name = "z-scores", 
        # row_dend_reorder = FALSE,
        row_dend_reorder = TRUE, 
        column_dend_reorder = FALSE, 
        show_column_dend = TRUE,
        show_heatmap_legend = FALSE,
        cluster_columns = TRUE,
        show_row_name = FALSE,
        column_names_rot = 45
        # width = ncol(tmt_lr_means_scaled)*unit(3, "cm"),
        # height = nrow(tmt_lr_means_scaled)*unit(3, "cm")
) 


heatmap_DEP_grouped_20 <- Heatmap(tmt_lr_means_scaled, name = "z-scores", 
             # row_dend_reorder = FALSE,
             row_dend_reorder = FALSE, 
             column_dend_reorder = FALSE, 
             show_column_dend = TRUE,
             show_heatmap_legend = FALSE,
             cluster_columns = TRUE,
             show_row_name = TRUE,
             row_split = 20
             # width = ncol(tmt_lr_means_scaled)*unit(3, "cm"),
             # height = nrow(tmt_lr_means_scaled)*unit(3, "cm")
)


HM <- draw(heatmap_DEP_grouped_20) 
r.dend <- row_dend(HM)
rcl.list <- row_order(HM)


list_DEP_by_cluster <- lapply(1:length(rcl.list), function(i){
        out <- data.frame(GeneID = rownames(tmt_lr_means_scaled[])[rcl.list[[i]]],
                          Cluster = paste0("cluster_", i),
                          stringsAsFactors = FALSE)
        return(out)
        }) %>%  
        do.call(rbind, .)
		