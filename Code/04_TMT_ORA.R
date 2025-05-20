# --------------------------
# ORA Analysis
# --------------------------

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







# UP regulated protens
## GO enrichment analysis for molecular functions
go_molec_lowIDHHG_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$IDHHG_LowGrade_up == 1],
                              'org.Hs.eg.db',
                              keyType = "UNIPROT",
                              ont = "MF",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              universe = dta_lowrest_sign$ID,
                              minGSSize = 5,
                              pool = TRUE)

go_molec_lowRTKPN_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKPN_LowGrade_up == 1],
                              'org.Hs.eg.db',
                              keyType = "UNIPROT",
                              ont = "MF",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              universe = dta_lowrest_sign$ID,
                              minGSSize = 5,
                              pool = TRUE)

go_molec_lowRTKCL_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKCL_LowGrade_up == 1],
                              'org.Hs.eg.db',
                              keyType = "UNIPROT",
                              ont = "MF",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              universe = dta_lowrest_sign$ID,
                              minGSSize = 5,
                              pool = TRUE)


go_molec_lowMES_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$MES_LowGrade_up == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "MF",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)



## Go enrichment analysis for biological processes
go_bio_lowIDHHG_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$IDHHG_LowGrade_up == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)


go_bio_lowRTKPN_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKPN_LowGrade_up == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)

go_bio_lowRTKCL_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKCL_LowGrade_up == 1],
                            'org.Hs.eg.db',
                            keyType = "UNIPROT",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = dta_lowrest_sign$ID,
                            minGSSize = 5,
                            pool = TRUE)


go_bio_lowMES_up <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$MES_LowGrade_up == 1],
                          'org.Hs.eg.db',
                          keyType = "UNIPROT",
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = dta_lowrest_sign$ID,
                          minGSSize = 5,
                          pool = TRUE)




# - 
# DOWN regulated proteins
## GO enrichment analysis for molecular functions
go_molec_lowIDHHG_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$IDHHG_LowGrade_down == 1],
                                 'org.Hs.eg.db',
                                 keyType = "UNIPROT",
                                 ont = "MF",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 universe = dta_lowrest_sign$ID,
                                 minGSSize = 5,
                                 pool = TRUE)

go_molec_lowRTKPN_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKPN_LowGrade_down == 1],
                                 'org.Hs.eg.db',
                                 keyType = "UNIPROT",
                                 ont = "MF",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 universe = dta_lowrest_sign$ID,
                                 minGSSize = 5,
                                 pool = TRUE)

go_molec_lowRTKCL_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKCL_LowGrade_down == 1],
                                 'org.Hs.eg.db',
                                 keyType = "UNIPROT",
                                 ont = "MF",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 universe = dta_lowrest_sign$ID,
                                 minGSSize = 5,
                                 pool = TRUE)


go_molec_lowMES_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$MES_LowGrade_down == 1],
                               'org.Hs.eg.db',
                               keyType = "UNIPROT",
                               ont = "MF",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               universe = dta_lowrest_sign$ID,
                               minGSSize = 5,
                               pool = TRUE)



## Go enrichment analysis for biological processes
go_bio_lowIDHHG_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$IDHHG_LowGrade_down == 1],
                               'org.Hs.eg.db',
                               keyType = "UNIPROT",
                               ont = "BP",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               universe = dta_lowrest_sign$ID,
                               minGSSize = 5,
                               pool = TRUE)


go_bio_lowRTKPN_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKPN_LowGrade_down == 1],
                               'org.Hs.eg.db',
                               keyType = "UNIPROT",
                               ont = "BP",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               universe = dta_lowrest_sign$ID,
                               minGSSize = 5,
                               pool = TRUE)

go_bio_lowRTKCL_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$RTKCL_LowGrade_down == 1],
                               'org.Hs.eg.db',
                               keyType = "UNIPROT",
                               ont = "BP",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               universe = dta_lowrest_sign$ID,
                               minGSSize = 5,
                               pool = TRUE)


go_bio_lowMES_down <- enrichGO(dta_lowrest_sign$ID[dta_lowrest_sign$MES_LowGrade_down == 1],
                             'org.Hs.eg.db',
                             keyType = "UNIPROT",
                             ont = "BP",
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH",
                             universe = dta_lowrest_sign$ID,
                             minGSSize = 5,
                             pool = TRUE)







## Variables 
network_layout      <- "kk"
n_categories_cnet   <- 6



########################
## Loop to replace UNIPORT-ID with GENE names 
rename_genes_ora <- function(x, df){
    
    if(!is.null(x)){
        
        temp <- slot(x, "result")[,"geneID"]
        
        for(i in df[["ID"]]) {     
            temp <- str_replace(temp, 
                                i, 
                                str_replace(df[df[["ID"]] == i, "name"], "_HUMAN", "")
            )
        }
        temp
    }
}


# -------------------------------------


# Low grade vs IDHHG
ORA_GO_bioProcesses_IDHHG_up                 <- go_bio_lowIDHHG_up
# ORA_GO_bioProcesses_IDHHG_up@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_IDHHG_up, dta_lowrest_sign)

ORA_GO_molecFunctions_IDHHG_up               <- go_molec_lowIDHHG_up
# ORA_GO_molecFunctions_IDHHG_up@result$geneID <- rename_genes_ora(ORA_GO_molecFunctions_IDHHG_up, dta_lowrest_sign)

# Low grade vs RTKPN
ORA_GO_bioProcesses_GBPN_up                 <- go_bio_lowRTKPN_up
ORA_GO_bioProcesses_GBPN_up@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_GBPN_up, dta_lowrest_sign)

ORA_GO_molecFunctions_GBPN_up               <- go_molec_lowRTKPN_up
ORA_GO_molecFunctions_GBPN_up@result$geneID <- rename_genes_ora(ORA_GO_molecFunctions_GBPN_up, dta_lowrest_sign)

# Low grade vs RTKCL
ORA_GO_bioProcesses_GBCL_up                 <- go_bio_lowRTKCL_up
ORA_GO_bioProcesses_GBCL_up@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_GBCL_up, dta_lowrest_sign)

ORA_GO_molecFunctions_GBCL_up               <- go_molec_lowRTKCL_up
ORA_GO_molecFunctions_GBCL_up@result$geneID <- rename_genes_ora(ORA_GO_molecFunctions_GBCL_up, dta_lowrest_sign)

# Low grade vs MES
ORA_GO_bioProcesses_GBMES_up                 <- go_bio_lowMES_up
ORA_GO_bioProcesses_GBMES_up@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_GBMES_up, dta_lowrest_sign)

ORA_GO_molecFunctions_GBMES_up               <- go_molec_lowMES_up
ORA_GO_molecFunctions_GBMES_up@result$geneID <- rename_genes_ora(ORA_GO_molecFunctions_GBMES_up, dta_lowrest_sign)





# Low grade vs IDHHG
ORA_GO_bioProcesses_IDHHG_down                 <- go_bio_lowIDHHG_down
ORA_GO_bioProcesses_IDHHG_down@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_IDHHG_down, dta_lowrest_sign)

ORA_GO_molecFunc_IDHHG_down               <- go_molec_lowIDHHG_down
# ORA_GO_molecFunc_IDHHG_down@result$geneID <- rename_genes_ora(ORA_GO_molecFunc_IDHHG_down, dta_lowrest_sign)


# Low grade vs RTKPN
ORA_GO_bioProcesses_GBPN_down                 <- go_bio_lowRTKPN_down
ORA_GO_bioProcesses_GBPN_down@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_GBPN_down, dta_lowrest_sign)

ORA_GO_molecFunc_GBPN_down               <- go_molec_lowRTKPN_down
ORA_GO_molecFunc_GBPN_down@result$geneID <- rename_genes_ora(ORA_GO_molecFunc_GBPN_down, dta_lowrest_sign)

# Low grade vs RTKCL
ORA_GO_bioProcesses_GBCL_down                 <- go_bio_lowRTKCL_down
ORA_GO_bioProcesses_GBCL_down@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_GBCL_down, dta_lowrest_sign)

ORA_GO_molecFunc_GBCL_down               <- go_molec_lowRTKCL_down
ORA_GO_molecFunc_GBCL_down@result$geneID <- rename_genes_ora(ORA_GO_molecFunc_GBCL_down, dta_lowrest_sign)

# Low grade vs MES
ORA_GO_bioProcesses_GBMES_down                 <- go_bio_lowMES_down
ORA_GO_bioProcesses_GBMES_down@result$geneID   <- rename_genes_ora(ORA_GO_bioProcesses_GBMES_down, dta_lowrest_sign)

ORA_GO_molecFunc_GBMES_down               <- go_molec_lowMES_down
ORA_GO_molecFunc_GBMES_down@result$geneID <- rename_genes_ora(ORA_GO_molecFunc_GBMES_down, dta_lowrest_sign)







