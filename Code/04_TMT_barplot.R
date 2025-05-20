
# -----------------------
# Barplot with number of DEP by contrasts
# -----------------------

## Create a temporary data.frame
dta_temp      <- as.data.frame(tmt_lr_of_interest)
## Create empty columns 
dta_temp[, contrasts_TMT_lr] <- NA
up_down <- data.frame(contrast = NULL, up = NULL, down = NULL)


## Looping through all contrast data set with the DEP
## 
for(i in list_of_tmt_lr)
{
    # row_name <- rownames(eval(parse(text = i)))
    dta_protein <- eval(parse(text = i))
    temp        <- gsub("proteins_", "" ,i)
    
    
    # Save the proteins that are up- and down-regulated separately
    dta_up   <- dta_protein[dta_protein[["logFC"]]>0,]
    dta_down <- dta_protein[dta_protein[["logFC"]]<0,]
    
    dta_temp[dta_temp$tmt_lr_of_interest %in% rownames(dta_up) , temp]   <- "Up"
    dta_temp[dta_temp$tmt_lr_of_interest %in% rownames(dta_down) , temp] <- "Down"
    
    ## Combine dataset 
    ## Count the number of up- & donw-regulated proteins by contrast
    up_down <- rbind(up_down, c(temp, 
                                sum(dta_temp[,temp] == "Up", na.rm = TRUE),
                                sum(dta_temp[,temp] == "Down", na.rm = TRUE)))
    
    
}


## Dataset with up and down regulated proteins by comparision
dta_up_down_lr <- dta_temp

## Change column names
colnames(up_down) <- c("contrast", "Up", "Down")
## Change orientation of data set from lonh to wide 
up_down <- up_down %>%                                  
    pivot_longer(c("Up", "Down"))

## Change dataset cosmetics
# Save values as numeric
up_down$value    <- as.numeric(up_down$value)
# Replace _ wiht / between contrats
up_down$contrast <- gsub("_", " / " ,up_down$contrast)
# Save order of contrast on X- Axis (otherwise it would be alphabetical)
up_down$contrast <- factor(up_down$contrast, levels = unique(up_down$contrast))
# Save order of Up and down - UP first otherwise it would be alphabetical)
up_down$name     <- factor(up_down$name,     levels = unique(up_down$name))

# Change order of the rows
up_down$contrast <- as.factor(up_down$contrast)
up_down$contrast <- factor(up_down$contrast, 
                           levels = c("IDHHG / LowGrade", "RTKPN / LowGrade", "RTKCL / LowGrade", "MES / LowGrade"))


up_down <- mutate(up_down, contrast = forcats::fct_recode(contrast, 
                                                          "IDH HG / Low Grade Glioma" = "IDHHG / LowGrade" ,
                                                          "GB PN / Low Grade Glioma" = "RTKPN / LowGrade",
                                                          "GB MES / Low Grade Glioma" =  "MES / LowGrade" ,
                                                          "GB CL / Low Grade Glioma" = "RTKCL / LowGrade"))

## Bar plot 
p_bar <- ggplot(up_down, aes(x = contrast, y = value, fill = name  )) + 
    geom_bar(stat = "identity", color = "black") + 
    theme_classic() +
    # scale_fill_brewer() +
    scale_fill_manual(values = c(col_up, col_down)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(0, 0, 0, 60),
          text = element_text(size = 22)) +
    ylab("Number of differentially expressed proteins") +
    xlab("") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(expand = c(0,0)) 
p_bar





