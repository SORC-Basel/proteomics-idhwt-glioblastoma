## -----------------------------------------
## GO Term identification
## -----------------------------------------
# This part is based on the cluster profiler package 

## New dataset with indicator if a protein was a significant DEP
dta_lowrest_sign         <- dta_tmt_low_rest_conf

## contrasts_TMT <- gsub("proteins_", "" ,list_of_tmt)
sign_tmt_lr      <- as.data.frame(tmt_lr_of_interest)

## Create empty columns 
sign_tmt_lr[,contrasts_TMT_lr] <- 0

## 
for(i in list_of_tmt_lr)
{
    row_name <- rownames(eval(parse(text = i)))
    # print(row_name)
    temp <- gsub("proteins_", "" ,i)
    
    sign_tmt_lr[sign_tmt_lr$tmt_lr_of_interest  %in% row_name , temp] <- 1
    
}

dta_lowrest_sign$IDHHG_LowGrade <- 0
dta_lowrest_sign$RTKPN_LowGrade <- 0
dta_lowrest_sign$RTKCL_LowGrade <- 0
dta_lowrest_sign$MES_LowGrade   <- 0

dta_lowrest_sign$IDHHG_LowGrade[dta_lowrest_sign$name %in% sign_tmt_lr$tmt_lr_of_interest[sign_tmt_lr$IDHHG_LowGrade == 1]] <- 1
dta_lowrest_sign$RTKPN_LowGrade[dta_lowrest_sign$name %in% sign_tmt_lr$tmt_lr_of_interest[sign_tmt_lr$RTKPN_LowGrade == 1]] <- 1
dta_lowrest_sign$RTKCL_LowGrade[dta_lowrest_sign$name %in% sign_tmt_lr$tmt_lr_of_interest[sign_tmt_lr$RTKCL_LowGrade == 1]] <- 1
dta_lowrest_sign$MES_LowGrade[dta_lowrest_sign$name %in% sign_tmt_lr$tmt_lr_of_interest[sign_tmt_lr$MES_LowGrade == 1]] <- 1



dta_lowrest_sign$IDHHG_LowGrade_up <- 0
dta_lowrest_sign$RTKPN_LowGrade_up <- 0
dta_lowrest_sign$RTKCL_LowGrade_up <- 0
dta_lowrest_sign$MES_LowGrade_up   <- 0

dta_lowrest_sign$IDHHG_LowGrade_up[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$IDHHG_LowGrade == "Up"]] <- 1
dta_lowrest_sign$RTKPN_LowGrade_up[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$RTKPN_LowGrade == "Up"]] <- 1
dta_lowrest_sign$RTKCL_LowGrade_up[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$RTKCL_LowGrade == "Up"]] <- 1
dta_lowrest_sign$MES_LowGrade_up[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$MES_LowGrade == "Up"]] <- 1



dta_lowrest_sign$IDHHG_LowGrade_down <- 0
dta_lowrest_sign$RTKPN_LowGrade_down <- 0
dta_lowrest_sign$RTKCL_LowGrade_down <- 0
dta_lowrest_sign$MES_LowGrade_down   <- 0

dta_lowrest_sign$IDHHG_LowGrade_down[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$IDHHG_LowGrade == "Down"]] <- 1
dta_lowrest_sign$RTKPN_LowGrade_down[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$RTKPN_LowGrade == "Down"]] <- 1
dta_lowrest_sign$RTKCL_LowGrade_down[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$RTKCL_LowGrade == "Down"]] <- 1
dta_lowrest_sign$MES_LowGrade_down[dta_lowrest_sign$name %in% dta_up_down_lr$tmt_lr_of_interest[dta_up_down_lr$MES_LowGrade == "Down"]] <- 1






## GO enrichment analysis for molecular functions
go_molec_lowIDHHG <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$IDHHG_LowGrade == 1],
                              'org.Hs.eg.db',
                              keyType = "UNIPROT",
                              ont = "MF",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              universe = dta_lowrest_sign$ID,
                              minGSSize = 5,
                              pool = TRUE)

go_molec_lowRTKPN <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKPN_LowGrade == 1],
                              'org.Hs.eg.db',
                              keyType = "UNIPROT",
                              ont = "MF",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              universe = dta_lowrest_sign$ID,
                              minGSSize = 5,
                              pool = TRUE)

go_molec_lowRTKCL <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKCL_LowGrade == 1],
                              'org.Hs.eg.db',
                              keyType = "UNIPROT",
                              ont = "MF",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              universe = dta_lowrest_sign$ID,
                              minGSSize = 5,
                              pool = TRUE)


go_molec_lowMES <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$MES_LowGrade == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "MF",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)



## Go enrichment analysis for biological processes
go_bio_lowIDHHG <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$IDHHG_LowGrade == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)


go_bio_lowRTKPN <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKPN_LowGrade == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)

go_bio_lowRTKCL <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKCL_LowGrade == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)


go_bio_lowMES <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$MES_LowGrade == 1],
                          'org.Hs.eg.db',
                          keyType = "UNIPROT",
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = dta_lowrest_sign$ID,
                          minGSSize = 5,
                          pool = TRUE)


# ## Go enrichment analysis for cellular components 
gsea_idhhg             <- data.frame("significant" = dta_lowrest_sign$IDHHG_LowGrade)
gsea_idhhg$logfc_names <- dta_lowrest_sign$ID
gsea_idhhg$logfc       <- fit_tmt_lr_2_full$coefficients[,"IDHHG - LowGrade"]
gsea_idhhg_vec         <- gsea_idhhg$logfc
names(gsea_idhhg_vec)  <- gsea_idhhg$logfc_names
gsea_idhhg_vec         <- sort(gsea_idhhg_vec, decreasing = TRUE)


gsea_rtkcl             <- data.frame("significant" = dta_lowrest_sign$RTKCL_LowGrade)
gsea_rtkcl$logfc_names <- dta_lowrest_sign$ID
gsea_rtkcl$logfc       <- fit_tmt_lr_2_full$coefficients[,"RTKCL - LowGrade"]
gsea_rtkcl_vec         <- gsea_rtkcl$logfc
names(gsea_rtkcl_vec)  <- gsea_rtkcl$logfc_names
gsea_rtkcl_vec         <- sort(gsea_rtkcl_vec, decreasing = TRUE)


gsea_rtkpn             <- data.frame("significant" = dta_lowrest_sign$RTKPN_LowGrade)
gsea_rtkpn$logfc_names <- dta_lowrest_sign$ID
gsea_rtkpn$logfc       <- fit_tmt_lr_2_full$coefficients[,"RTKPN - LowGrade"]
gsea_rtkpn_vec         <- gsea_rtkpn$logfc
names(gsea_rtkpn_vec)  <- gsea_rtkpn$logfc_names
gsea_rtkpn_vec         <- sort(gsea_rtkpn_vec, decreasing = TRUE)

gsea_mes             <- data.frame("significant" = dta_lowrest_sign$MES_LowGrade)
gsea_mes$logfc_names <- dta_lowrest_sign$ID
gsea_mes$logfc       <- fit_tmt_lr_2_full$coefficients[,"MES - LowGrade"]
gsea_mes_vec         <- gsea_mes$logfc
names(gsea_mes_vec)  <- gsea_mes$logfc_names
gsea_mes_vec         <- sort(gsea_mes_vec, decreasing = TRUE)





gsea_idhhg_p             <- data.frame("significant" = dta_lowrest_sign$IDHHG_LowGrade)
gsea_idhhg_p$logfc_names <- dta_lowrest_sign$ID
gsea_idhhg_p$logfc       <- fit_tmt_lr_2_full$coefficients[,"IDHHG - LowGrade"]
gsea_idhhg_p$p_value     <- fit_tmt_lr_2_full$p.value[,"IDHHG - LowGrade"]
gsea_idhhg_p_vec         <- sign(gsea_idhhg_p$logfc) * -1 * log10(gsea_idhhg_p$p_value)
names(gsea_idhhg_p_vec)  <- gsea_idhhg_p$logfc_names
gsea_idhhg_p_vec         <- sort(gsea_idhhg_p_vec, decreasing = TRUE)


gsea_rtkcl_p             <- data.frame("significant" = dta_lowrest_sign$RTKCL_LowGrade)
gsea_rtkcl_p$logfc_names <- dta_lowrest_sign$ID
gsea_rtkcl_p$logfc       <- fit_tmt_lr_2_full$coefficients[,"RTKCL - LowGrade"]
gsea_rtkcl_p$p_value     <- fit_tmt_lr_2_full$p.value[,"RTKCL - LowGrade"]
gsea_rtkcl_p_vec         <- sign(gsea_rtkcl_p$logfc) * -1 * log10(gsea_rtkcl_p$p_value)
names(gsea_rtkcl_p_vec)  <- gsea_rtkcl_p$logfc_names
gsea_rtkcl_p_vec         <- sort(gsea_rtkcl_p_vec, decreasing = TRUE)


gsea_rtkpn_p             <- data.frame("significant" = dta_lowrest_sign$RTKPN_LowGrade)
gsea_rtkpn_p$logfc_names <- dta_lowrest_sign$ID
gsea_rtkpn_p$logfc       <- fit_tmt_lr_2_full$coefficients[,"RTKPN - LowGrade"]
gsea_rtkpn_p$p_value     <- fit_tmt_lr_2_full$p.value[,"RTKPN - LowGrade"]
gsea_rtkpn_p_vec         <- sign(gsea_rtkpn_p$logfc) * -1 * log10(gsea_rtkpn_p$p_value)
names(gsea_rtkpn_p_vec)  <- gsea_rtkpn_p$logfc_names
gsea_rtkpn_p_vec         <- sort(gsea_rtkpn_p_vec, decreasing = TRUE)


gsea_mes_p             <- data.frame("significant" = dta_lowrest_sign$MES_LowGrade)
gsea_mes_p$logfc_names <- dta_lowrest_sign$ID
gsea_mes_p$logfc       <- fit_tmt_lr_2_full$coefficients[,"MES - LowGrade"]
gsea_mes_p$p_value     <- fit_tmt_lr_2_full$p.value[,"MES - LowGrade"]
gsea_mes_p_vec         <- sign(gsea_mes_p$logfc) * -1 * log10(gsea_mes_p$p_value)
names(gsea_mes_p_vec)  <- gsea_mes_p$logfc_names
gsea_mes_p_vec         <- sort(gsea_mes_p_vec, decreasing = TRUE)






# GO - GSEA on signed P-Values
go_gsea <- function(x, ont_type = "BP"){
    gseGO(x,
          'org.Hs.eg.db',
          keyType = "UNIPROT",
          ont = ont_type,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          maxGSSize = 3000,
          eps = 1e-30
          # universe = dta_lowrest_sign$ID,
          # minGSSize = 5,
          # pool = TRUE
    )
    
}




## GSEA - Go enrichment analysis for biological processes
gsea_bio_lowIDHHG <- go_gsea(gsea_idhhg_vec, ont_type = "BP")
gsea_bio_lowRTKCL <- go_gsea(gsea_rtkcl_vec, ont_type = "BP")
gsea_bio_lowRTKPN <- go_gsea(gsea_rtkpn_vec, ont_type = "BP")
gsea_bio_lowMES <- go_gsea(gsea_mes_vec, ont_type = "BP")

## GSEA - Go enrichment analysis for molecular functions
gsea_molec_lowIDHHG <- go_gsea(gsea_idhhg_vec, ont_type = "MF")
gsea_molec_lowRTKCL <- go_gsea(gsea_rtkcl_vec, ont_type = "MF")
gsea_molec_lowRTKPN <- go_gsea(gsea_rtkpn_vec, ont_type = "MF")
gsea_molec_lowMES <- go_gsea(gsea_mes_vec, ont_type = "MF")

## GSEA - Go enrichment analysis for biological processes
gsea_bio_lowIDHHG_p <- go_gsea(gsea_idhhg_p_vec, ont_type = "BP")
gsea_bio_lowRTKCL_p <-  go_gsea(gsea_rtkcl_p_vec, ont_type = "BP")
gsea_bio_lowRTKPN_p <- go_gsea(gsea_rtkpn_p_vec, ont_type = "BP"
gsea_bio_lowMES_p <- go_gsea(gsea_mes_p_vec, ont_type = "BP")

## GSEA - Go enrichment analysis for molecular functions
gsea_molec_lowIDHHG_p <- go_gsea(gsea_idhhg_p_vec, ont_type = "MF")
gsea_molec_lowRTKCL_p <- go_gsea(gsea_rtkcl_p_vec, ont_type = "MF")
gsea_molec_lowRTKPN_p <- go_gsea(gsea_rtkpn_p_vec, ont_type = "MF")
gsea_molec_lowMES_p <- go_gsea(gsea_mes_p_vec, ont_type = "MF")







# In both LogFC and P-Value
go_bio_gsea_idhhg_p_both <- gsea_bio_lowIDHHG_p@result[gsea_bio_lowIDHHG_p@result$ID %in% gsea_bio_lowIDHHG@result$ID,]
go_bio_gsea_idhhg_both   <- gsea_bio_lowIDHHG@result[gsea_bio_lowIDHHG@result$ID %in% gsea_bio_lowIDHHG_p@result$ID,]

go_bio_gsea_rtkcl_p_both <- gsea_bio_lowRTKCL_p@result[gsea_bio_lowRTKCL_p@result$ID %in% gsea_bio_lowRTKCL@result$ID,]
go_bio_gsea_rtkcl_both   <- gsea_bio_lowRTKCL@result[gsea_bio_lowRTKCL@result$ID %in% gsea_bio_lowRTKCL_p@result$ID,]

go_bio_gsea_rtkpn_p_both <- gsea_bio_lowRTKPN_p@result[gsea_bio_lowRTKPN_p@result$ID %in% gsea_bio_lowRTKPN@result$ID,]
go_bio_gsea_rtkpn_both   <- gsea_bio_lowRTKPN@result[gsea_bio_lowRTKPN@result$ID %in% gsea_bio_lowRTKPN_p@result$ID,]

go_bio_gsea_mes_p_both   <- gsea_bio_lowMES_p@result[gsea_bio_lowMES_p@result$ID %in% gsea_bio_lowMES@result$ID,]
go_bio_gsea_mes_both     <- gsea_bio_lowMES@result[gsea_bio_lowMES@result$ID %in% gsea_bio_lowMES_p@result$ID,]


col_names <- c("ID", "Description", 
               "setSize_Pvalue", "enrichmentScore_Pvalue", "NES_pvalue", "pvalue_Pvalue", "p.adjust_Pvalue", 
               "qvalue_Pvalue", "rank_Pvalue", "leading_edge_Pvalue", "core_enrichment_Pvalue",
               "setSize_logfc", "enrichmentScore_logfc", "NES_pvalue", "pvalue_logfc", "p.adjust_logfc", 
               "qvalue_logfc", "rank_logfc", "leading_edge_logfc", "core_enrichment_logfc"
)

go_bio_gsea_idhh_merged  <- merge(go_bio_gsea_idhhg_p_both, go_bio_gsea_idhhg_both, by = c("ID", "Description"))
go_bio_gsea_rtkcl_merged <- merge(go_bio_gsea_rtkcl_p_both, go_bio_gsea_rtkcl_both, by = c("ID", "Description"))
go_bio_gsea_rtkpn_merged <- merge(go_bio_gsea_rtkpn_p_both, go_bio_gsea_rtkpn_both, by = c("ID", "Description"))
go_bio_gsea_mes_merged   <- merge(go_bio_gsea_mes_p_both, go_bio_gsea_mes_both, by = c("ID", "Description"))

colnames(go_bio_gsea_idhh_merged) <- col_names
colnames(go_bio_gsea_rtkcl_merged) <- col_names
colnames(go_bio_gsea_rtkpn_merged) <- col_names
colnames(go_bio_gsea_mes_merged) <- col_names


go_molec_gsea_idhhg_p_both <- gsea_molec_lowIDHHG_p@result[gsea_molec_lowIDHHG_p@result$ID %in% gsea_molec_lowIDHHG@result$ID,]
go_molec_gsea_idhhg_both   <- gsea_molec_lowIDHHG@result[gsea_molec_lowIDHHG@result$ID %in% gsea_molec_lowIDHHG_p@result$ID,]

go_molec_gsea_rtkcl_p_both <- gsea_molec_lowRTKCL_p@result[gsea_molec_lowRTKCL_p@result$ID %in% gsea_molec_lowRTKCL@result$ID,]
go_molec_gsea_rtkcl_both   <- gsea_molec_lowRTKCL@result[gsea_molec_lowRTKCL@result$ID %in% gsea_molec_lowRTKCL_p@result$ID,]

go_molec_gsea_rtkpn_p_both <- gsea_molec_lowRTKPN_p@result[gsea_molec_lowRTKPN_p@result$ID %in% gsea_molec_lowRTKPN@result$ID,]
go_molec_gsea_rtkpn_both   <- gsea_molec_lowRTKPN@result[gsea_molec_lowRTKPN@result$ID %in% gsea_molec_lowRTKPN_p@result$ID,]

go_molec_gsea_mes_p_both   <- gsea_molec_lowMES_p@result[gsea_molec_lowMES_p@result$ID %in% gsea_molec_lowMES@result$ID,]
go_molec_gsea_mes_both     <- gsea_molec_lowMES@result[gsea_molec_lowMES@result$ID %in% gsea_molec_lowMES_p@result$ID,]


go_molec_gsea_idhh_merged  <- merge(go_molec_gsea_idhhg_p_both, go_molec_gsea_idhhg_both, by = c("ID", "Description"))
go_molec_gsea_rtkcl_merged <- merge(go_molec_gsea_rtkcl_p_both, go_molec_gsea_rtkcl_both, by = c("ID", "Description"))
go_molec_gsea_rtkpn_merged <- merge(go_molec_gsea_rtkpn_p_both, go_molec_gsea_rtkpn_both, by = c("ID", "Description"))
go_molec_gsea_mes_merged   <- merge(go_molec_gsea_mes_p_both, go_molec_gsea_mes_both, by = c("ID", "Description"))

colnames(go_molec_gsea_idhh_merged) <- col_names
colnames(go_molec_gsea_rtkcl_merged) <- col_names
colnames(go_molec_gsea_rtkpn_merged) <- col_names
colnames(go_molec_gsea_mes_merged) <- col_names









## Variables 
network_layout      <- "kk"
n_categories_cnet   <- 6

## Simplify GO terms
# Enrichment
simple_go_bio_lowIDHHG     <- simplify(go_bio_lowIDHHG)
# simple_go_molec_lowIDHHG   <- simplify(go_molec_lowIDHHG)

simple_go_bio_lowRTKPN     <- simplify(go_bio_lowRTKPN)
simple_go_molec_lowRTKPN   <- simplify(go_molec_lowRTKPN)

simple_go_bio_lowRTKCL     <- simplify(go_bio_lowRTKCL)
simple_go_molec_lowRTKCL   <- simplify(go_molec_lowRTKCL)

simple_go_bio_lowMES       <- simplify(go_bio_lowMES)
simple_go_molec_lowMES     <- simplify(go_molec_lowMES)



# GSEA - Log FC
simple_gsea_bio_lowIDHHG     <- simplify(gsea_bio_lowIDHHG)
simple_gsea_molec_lowIDHHG   <- simplify(gsea_molec_lowIDHHG)

simple_gsea_bio_lowRTKPN     <- simplify(gsea_bio_lowRTKPN)
simple_gsea_molec_lowRTKPN   <- simplify(gsea_molec_lowRTKPN)

simple_gsea_bio_lowRTKCL     <- simplify(gsea_bio_lowRTKCL)
simple_gsea_molec_lowRTKCL   <- simplify(gsea_molec_lowRTKCL)

simple_gsea_bio_lowMES       <- simplify(gsea_bio_lowMES)
simple_gsea_molec_lowMES     <- simplify(gsea_molec_lowMES)




########################
## Loop to replace UNIPORT-ID with GENE names 
rename_genes <- function(x, df){
    temp <- slot(x, "result")[,"geneID"]
    
    for(i in df[["ID"]]) {     
        temp <- str_replace(temp, 
                            i, 
                            str_replace(df[df[["ID"]] == i, "name"], "_HUMAN", "")
        )
    }
    temp
}



## Create gene list with foldchange
foldchange_lr        <- fit_tmt_lr_2$coefficients
names(foldchange_lr) <- str_replace(rownames(foldchange_lr), "_HUMAN", "") 



## Function to creratre GSEA cnetplots
go_cnetplot <- function(x, by.count = TRUE){
    
    # if(nrow(slot(x, "result")) == 0)
    if(nrow(slot(x, "result")) == 0) {
        print("No enriched path identified")
    } else {
        
        # Low grade vs IDHHG
        go_gene                        <- x
        go_gene@result$geneID <- rename_genes(x, dta_lowrest_sign)
        
        # Selection of categories 
        go_results <- slot(x, "result")
        # chosen by enrichment score
        if(by.count){
            go_results <- go_results[order(go_results$Count, decreasing = TRUE),]
            go_results <- go_results[go_results$p.adjust<0.05,]
        }
        
        # List of shown categories
        go_list_cnet <- go_results[1:n_categories_cnet, "Description"]
        
        
        print(go_list_cnet)
        p <- cnetplot(go_gene, 
                      showCategory       = go_list_cnet,
                      cex_label_category = 0.8, 
                      cex_label_gene     = 0.4,
                      colorEdge          = TRUE,
                      color_category     = 'firebrick',
                      color_gene         = 'steelblue',
                      layout             = network_layout,
                      foldChange         = foldchange_lr
        ) + 
            guides(edge_colour = "none") +
            guides(size = "none") +
            scale_colour_gradient2(name = "fold change", low = "blue",
                                   mid = "white", high = "red",
                                   guide = guide_colorbar(order = 2),
                                   midpoint = 0,
                                   limits = c(min(foldchange_lr), max(foldchange_lr)),
                                   breaks = scales::pretty_breaks(n=5)(min(foldchange_lr):max(foldchange_lr)))
        
        
        p
    }
}




p_lowIDHHG_bio   <- go_cnetplot(go_bio_lowIDHHG)
# p_lowIDHHG_molec <- go_cnetplot(go_molec_lowIDHHG)

p_lowRTKPN_bio   <- go_cnetplot(go_bio_lowRTKPN)
# p_lowRTKPN_molec <- go_cnetplot(go_molec_lowRTKPN)

# p_lowRTKCL_bio   <- go_cnetplot(go_bio_lowRTKCL)
# p_lowRTKCL_molec <- go_cnetplot(go_molec_lowRTKCL)

# p_lowMES_bio     <- go_cnetplot(go_bio_lowMES)
# p_lowMES_molec   <- go_cnetplot(go_molec_lowMES)








########################
## Loop to replace UNIPORT-ID with GENE names 
rename_genes_gsea <- function(x, df){
    temp <- slot(x, "result")[,"core_enrichment"]
    
    for(i in df[["ID"]]) {     
        temp <- str_replace(temp, 
                            i, 
                            str_replace(df[df[["ID"]] == i, "name"], "_HUMAN", "")
        )
    }
    temp
}


## Function to creratre GSEA cnetplots
gsea_cnetplot <- function(x, by.enrichment = TRUE){
    
    # if(nrow(slot(x, "result")) == 0)
    if(nrow(slot(x, "result")) == 0) {
        print("No enriched path identified")
    } else {
        
        # Low grade vs IDHHG
        gsea_gene                        <- x
        gsea_gene@result$core_enrichment <- rename_genes_gsea(x, dta_lowrest_sign)
        
        # Selection of categories 
        gsea_results <- slot(x, "result")
        # chosen by enrichment score
        if(by.enrichment){
            gsea_results <- gsea_results[order(abs(gsea_results$enrichmentScore), decreasing = TRUE),]
            gsea_results <- gsea_results[gsea_results$p.adjust<0.05,]
        }
        
        # List of shown categories
        gsea_list_cnet <- gsea_results[1:n_categories_cnet, "Description"]
        
        
        print(gsea_list_cnet)
        p <- cnetplot(gsea_gene, 
                      showCategory       = gsea_list_cnet,
                      cex_label_category = 0.8, 
                      cex_label_gene     = 0.4,
                      colorEdge          = TRUE,
                      color_category     = 'firebrick',
                      color_gene         = 'steelblue',
                      layout             = network_layout,
                      foldChange         = foldchange_lr
        ) + 
            guides(edge_colour = "none") +
            guides(size = "none") +
            scale_colour_gradient2(name = "fold change", low = "blue",
                                   mid = "white", high = "red",
                                   guide = guide_colorbar(order = 2),
                                   midpoint = 0,
                                   limits = c(min(foldchange_lr), max(foldchange_lr)),
                                   breaks = scales::pretty_breaks(n=5)(min(foldchange_lr):max(foldchange_lr)))
        
        
        p
    }
}




cnet_gsea_lowIDHHG_bio   <- gsea_cnetplot(gsea_bio_lowIDHHG, by.enrichment = FALSE)
cnet_gsea_lowIDHHG_molec <- gsea_cnetplot(gsea_molec_lowIDHHG, by.enrichment = FALSE)

cnet_gsea_lowRTKPN_bio   <- gsea_cnetplot(gsea_bio_lowRTKPN, by.enrichment = FALSE)
cnet_gsea_lowRTKPN_molec <- gsea_cnetplot(gsea_molec_lowRTKPN, by.enrichment = FALSE)

cnet_gsea_lowRTKCL_bio   <- gsea_cnetplot(gsea_bio_lowRTKCL, by.enrichment = FALSE)
cnet_gsea_lowRTKCL_molec <- gsea_cnetplot(gsea_molec_lowRTKCL, by.enrichment = FALSE)

cnet_gsea_lowMES_bio     <- gsea_cnetplot(gsea_bio_lowMES, by.enrichment = FALSE)
cnet_gsea_lowMES_molec   <- gsea_cnetplot(gsea_molec_lowMES, by.enrichment = FALSE)



# Cnetplot - P-Values
cnet_gsea_lowIDHHG_bio_p   <- gsea_cnetplot(gsea_bio_lowIDHHG_p, by.enrichment = FALSE)
cnet_gsea_lowIDHHG_molec_p <- gsea_cnetplot(gsea_molec_lowIDHHG_p, by.enrichment = FALSE)

cnet_gsea_lowRTKPN_bio_p   <- gsea_cnetplot(gsea_bio_lowRTKPN_p, by.enrichment = FALSE)
cnet_gsea_lowRTKPN_molec_p <- gsea_cnetplot(gsea_molec_lowRTKPN_p, by.enrichment = FALSE)

cnet_gsea_lowRTKCL_bio_p   <- gsea_cnetplot(gsea_bio_lowRTKCL_p, by.enrichment = FALSE)
cnet_gsea_lowRTKCL_molec_p <- gsea_cnetplot(gsea_molec_lowRTKCL_p, by.enrichment = FALSE)

cnet_gsea_lowMES_bio_p     <- gsea_cnetplot(gsea_bio_lowMES_p, by.enrichment = FALSE)
cnet_gsea_lowMES_molec_p   <- gsea_cnetplot(gsea_molec_lowMES_p, by.enrichment = FALSE)



figure_width  <- 14
figure_height <- 14






cnet_gsea_lowIDHHG_bio + ggtitle("IDH HG - Low grade glioma - Biological processes")  + 
    cnet_gsea_lowRTKPN_bio + ggtitle("GB PN - Low grade glioma - Biological processes")  +
    cnet_gsea_lowRTKCL_bio + ggtitle("GB CL - Low grade glioma - Biological processes") +
    cnet_gsea_lowMES_bio + ggtitle("GB MES - Low grade glioma - Biological processes") 


cnet_gsea_lowRTKCL_molec + ggtitle("RTKCL - Low grade glioma - Molecular functions")  + 
    cnet_gsea_lowMES_molec + ggtitle("MES - Low grade glioma - Molecular functions")  +
    cnet_gsea_lowMES_molec + ggtitle("MES - Low grade glioma - Molecular functions") +
    cnet_gsea_lowMES_molec + ggtitle("MES - Low grade glioma - Molecular functions") 


