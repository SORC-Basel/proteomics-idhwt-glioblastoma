



## ---------------------------------------------------------
## Vulcano plots for TMT proteomics
## ---------------------------------------------------------

# Function for ggplot vulcano plot - TMT
ggplot_vulcano <- function(model, coefficient, p_value = c("unadjusted", "adjusted"), labels = FALSE ){
    
    # Vulcanoplot - Unadjusted p-value
    limma.results <- topTable(model, coef = coefficient, n = Inf)
    
    limma.results$log.P.Val  = -log10(limma.results$P.Value)
    limma.results$prot.names <- str_replace(rownames(limma.results), "_HUMAN", "")
    limma.results$prot.names <- str_replace(limma.results$prot.names, "VIME_CON-HUMAN,VIME", "VIME")
    limma.results$prot.names <- str_replace(limma.results$prot.names, "LV746;LV743_HUMAN", "IGLV746")
    
    
    print_name <- switch(coefficient, 
                    "IDHHG - LowGrade"  = "IDH HG - Low Grade Glioma",
                    "RTKPN - LowGrade"  = "GB PN - Low Grade Glioma",
                    "RTKCL - LowGrade"  = "GB CL - Low Grade Glioma",
                    "MES - LowGrade"    = "GB MES - Low Grade Glioma" )
    

    
    if(p_value == "unadjusted"){
        limma.results.sign   =  limma.results[limma.results$adj.P.Val < 0.05,]
        
        p <- ggplot(limma.results, aes(x = logFC, y =log.P.Val )) +
            geom_point(size=0.8 )+
            #ggtitle(paste("TMT proteomics:", coefficient)) + 
            ggtitle(paste(print_name)) + 
            theme_bw(base_size = 16) + # change theme
            xlab(expression("log2 Fold Change")) + # x-axis label
            ylab(expression(" -log10(P-value)")) + # y-axis label
            
            # geom_vline(xintercept = c(-1,1), colour = "blue") + # Add fold change cutoffs
            # geom_hline(yintercept = -log10(0.05), colour = "red") + # Add significance cutoffs
            geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
            theme(legend.position = "none")
        
        # scale_colour_gradient(low = "black", high = "black", guide = FALSE)
        # geom_text_repel(data=subset(limma.results, abs(logFC)>1 & log.adj.P.Val > 1.5),
        #                 aes( logFC, log.adj.P.Val ,label=prot.names)) # add gene label
        # if(labels){  p <- p + geom_text_repel(data=subset(limma.results, abs(logFC)>0.75 & adj.P.Val < 0.05), aes( logFC, log.P.Val ,label=prot.names)) }
        if(labels) {
            # Subset the data for significant results with p < 0.05
            significant_results <- subset(limma.results, adj.P.Val < 0.05)
            
            # Sort by logFC, and handle the case where there are fewer than 3 entries in either group
            top_3 <- head(significant_results[order(-significant_results$logFC), ], min(3, nrow(significant_results)))
            bottom_3 <- head(significant_results[order(significant_results$logFC), ], min(3, nrow(significant_results)))
            
            # Combine top and bottom groups, ensuring that they only include available data
            labeled_results <- unique(rbind(top_3, bottom_3))
            
            # Add the labels only for these top and bottom entries
            p <- p + geom_text_repel(data=labeled_results, aes(logFC, log.P.Val, label=prot.names))
        }
        
        if(nrow(limma.results.sign) > 0){
            p <- p + geom_point(data = limma.results.sign[limma.results.sign["logFC"] > 0,],
                                aes( x = logFC, y = log.P.Val), color = col_up,
                                size = 1.5)
            
            p <- p + geom_point(data = limma.results.sign[limma.results.sign["logFC"] < 0,],
                                aes( x = logFC, y = log.P.Val), color = col_down,
                                size = 1.5)
        }
    }
    
    if(p_value == "adjusted"){
        
        limma.results$log.P.Val  = -log10(limma.results$adj.P.Val)
        
        p <- ggplot(limma.results, aes(x = logFC, y = log.P.Val )) +
            geom_point(size=0.8 )+
            ggtitle(paste("TMT:", coefficient)) + 
            theme_bw(base_size = 16) + # change theme
            xlab(expression("log2 Fold Change")) + # x-axis label
            ylab(expression(" -log10(adjust. P-value)")) + # y-axis label
            
            # geom_vline(xintercept = c(-1,1), colour = "blue") + # Add fold change cutoffs
            geom_hline(yintercept = -log10(0.05), colour = "red") + # Add significance cutoffs
            geom_vline(xintercept = 0, colour = "black") # Add 0 lines
        
        # scale_colour_gradient(low = "black", high = "black", guide = FALSE)
        # geom_text_repel(data=subset(limma.results, abs(logFC)>1 & log.adj.P.Val > 1.5),
        #                 aes( logFC, log.adj.P.Val ,label=prot.names)) # add gene label
        if(labels){  p <- p + geom_text_repel(data=subset(limma.results, abs(logFC)>0.75 & adj.P.Val < 0.05), aes( logFC, log.P.Val, label = prot.names)) }
        
        
    }
    
    
    p <- p + scale_x_continuous(breaks = seq(-2,2,0.5))
    p
    
    # ggsave(paste0(path_directory,"Figures/vulcanoplot_tmt_", gsub(" - ","_",coefficient), ".png")
    #   ,p)
    
}


## ---------------------------------------------------------
## GO analysis function  -- uses protti package 
## ---------------------------------------------------------
go_analysis <- function(contrast = NULL, type = c("molecular", "biological", "cellular")){
  
  # go_type <- switch(type,
  #                   "molecular"  = "go_molecular_function",
  #                   "biological" = "go_biological_process", 
  #                   "cellular"   = "go_cellular_compartment",
  #                   "go_molecular_function")
  # print(go_type)
  
  dta_contrast <- sign_tmt[,c("tmt_of_interest_dt",contrast)]
  dta_contrast <- dta_contrast[dta_contrast[[contrast]] == 1,]
  
  if(nrow(dta_contrast) == 0) stop("No singificant proteins for this contrats")
  
  dta_tmt_sign <- dta_full_tmt
  dta_tmt_sign$is_sign <- FALSE
  dta_tmt_sign$is_sign[dta_tmt_sign$name %in% dta_contrast$tmt_of_interest_dt] <- 1 
  
  
  
  dta_tmt_sign <- merge(dta_tmt_sign , 
                        uniprot_list[,c("id", 
                                        "go_molecular_function", 
                                        "go_biological_process",
                                        "go_cellular_compartment" )],
                        by.x = "ID", 
                        by.y = "id")
  
  
  
  if(type == "molecular"){
    print("Go aspect: molecular function")
    
    go_output <- go_enrichment(
      dta_tmt_sign,
      protein_id = ID,
      is_significant = is_sign,
      go_annotations_uniprot = go_molecular_function,
      plot_cutoff = "adj_pval top10")
    
  } else if (type == "biological"){
    print("Go aspect: biological process")
    
    go_output <- go_enrichment(
      dta_tmt_sign,
      protein_id = ID,
      is_significant = is_sign,
      go_annotations_uniprot = go_biological_process)
    
  } else if (type == "cellular"){
    print("Go aspect: cellular component")
    
    go_output <- go_enrichment(
      dta_tmt_sign,
      protein_id = ID,
      is_significant = is_sign,
      go_annotations_uniprot = go_cellular_compartment)
    
  }
  
  go_output
  
}






## ---------------------------------------------------------
## GO analysis function  -- uses clusterprofiler package 
## ---------------------------------------------------------
go_cluster <- function(contrast = NULL, go_type = c("molecular", "biological", "cellular")){
  
  
  dta_contrast <- sign_tmt[,c("tmt_of_interest_dt",contrast)]
  dta_contrast <- dta_contrast[dta_contrast[[contrast]] == 1,]
  
  if(nrow(dta_contrast) == 0) stop("No singificant proteins for this contrats")
  
  dta_tmt_sign <- dta_full_tmt
  dta_tmt_sign$is_sign <- FALSE
  dta_tmt_sign$is_sign[dta_tmt_sign$name %in% dta_contrast$tmt_of_interest_dt] <- 1 
  
  
  go_ont <- switch(go_type,
                   "molecular"  = "MF",
                   "biological" = "BP", 
                   "cellular"   = "CC",
                   "MF"
  )
  
  print(paste("Go aspect:", go_type, "function"))
  
  go_output <- enrichGO(dta_tmt_sign$ID[dta_tmt_sign$is_sign == 1],
                        'org.Hs.eg.db',
                        keyType = "UNIPROT",
                        ont = go_ont,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        universe = dta_tmt_sign$ID,
                        minGSSize = 5,
                        pool = TRUE
  )
  
  go_output
}





## ---------------------------------------------------------
## Dotplot function  -- uses clusterprofiler package 
## ---------------------------------------------------------
dot_plot_cluster <- function(cluster_molecular, cluster_biological, cluster_cellular){
  
  dot_molec <- as.data.frame(cluster_molecular)
  dot_molec <- dot_molec %>% 
    mutate(type = "Molecular")  %>% 
    slice_min(p.adjust, n=10,  with_ties = F)
  
  dot_bio <- as.data.frame(cluster_biological)
  dot_bio <- dot_bio %>% 
    mutate(type = "Biological")  %>% 
    slice_min(p.adjust, n=10,  with_ties = F)
  
  dot_cell <- as.data.frame(cluster_cellular)
  dot_cell <- dot_cell %>% 
    mutate(type = "Cellular")  %>% 
    slice_min(p.adjust, n=10,  with_ties = F)
  
  dot_df <- rbind(dot_molec, dot_bio, dot_cell)
  
  for(i in 1:nrow(dot_df)){
    dot_df$GeneRatio[i] <- eval(parse(text = dot_df$GeneRatio[i]))  
  }
  
  dot_df$GeneRatio <- as.numeric(dot_df$GeneRatio)
  
  
  contrast_string <- str_split(deparse(substitute(cluster_molecular)), "_", simplify = TRUE)
  
  contrast_title <- paste("Gene Ontology functional analysis -", contrast_string[3], "/" ,contrast_string[4] )
  
  
    ggplot(dot_df, aes(x = GeneRatio, y = forcats::fct_reorder(Description, GeneRatio))) + 
    geom_point(aes(color = p.adjust, size = Count)) +
    theme_bw(base_size = 14) +
    # scale_colour_gradient(limits=c(0, plyr::round_any(max(dot_df$p.adjust), 0.001, f = ceiling)), 
    #                       low="red", high = "blue") +
    scale_colour_gradient(limits=c(0, 0.05), low="red", high = "blue") +
    ylab(NULL) +
    labs(color = "Adjusted P-value") +
    labs(size  = "Proteins") +
    ggtitle(contrast_title) +
    # scale_x_continuous(limits = c(0, 0.12), breaks = seq(0,0.12,0.01), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, plyr::round_any(max(dot_df$GeneRatio), 0.01, f = ceiling)), 
                       breaks = seq(0, plyr::round_any(max(dot_df$GeneRatio), 0.01, f = ceiling),0.01), 
                       expand = c(0,0)) +
    facet_grid(vars(type), scales = "free")
  
  
  
}