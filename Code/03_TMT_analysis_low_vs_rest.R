# -----------------------
# TMT - Analysis
# -----------------------
## Create design matrix
dta_tmt_low_rest_conf <- dta_tmt_log

# Contrats Low vs High
cond_tmt_low_rest   = as.factor(c(rep("LowGrade",length(cols_IDHLG)),
                                  rep("IDHHG", length(cols_IDHHG)),
                                  rep("LowGrade",  length(cols_IDHTERT)),
                                  rep("LowGrade",  length(cols_IDHCODEL)),
                                  rep("RTKPN", length(cols_RTKPN)),
                                  rep("RTKCL", length(cols_RTKCL)),
                                  rep("MES", length(cols_MES)) ))


## Desgin matrices 
# No confounders
design_tmt_low_rest = model.matrix(~ 0 + cond_tmt_low_rest) 

# Univariable
design_tmt_low_rest_age       = model.matrix(~ 0 + cond_tmt_low_rest + age )
design_tmt_low_rest_sex       = model.matrix(~ 0 + cond_tmt_low_rest + sex )
design_tmt_low_rest_steroids  = model.matrix(~ 0 + cond_tmt_low_rest + steroids )
design_tmt_low_rest_therapy   = model.matrix(~ 0 + cond_tmt_low_rest + therapy_any )
design_tmt_low_rest_epilepsy  = model.matrix(~ 0 + cond_tmt_low_rest + epilepsy )
design_tmt_low_rest_hg2a      = model.matrix(~ 0 + cond_tmt_low_rest + hg2a )
design_tmt_low_rest_sept3     = model.matrix(~ 0 + cond_tmt_low_rest + sept3 )
design_tmt_low_rest_ptprc     = model.matrix(~ 0 + cond_tmt_low_rest + ptprc )

# Base Model
design_tmt_low_rest_base = model.matrix(~ 0 + cond_tmt_low_rest + age + sex + steroids + therapy_any + 
                                            epilepsy)


# Full Model
design_tmt_low_rest_full = model.matrix(~ 0 + cond_tmt_low_rest + age + sex + steroids + therapy_any + 
                                              epilepsy + hg2a + sept3 + ptprc )

colnames(design_tmt_low_rest) = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest))



# Univariable
colnames(design_tmt_low_rest_age)      = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_age))
colnames(design_tmt_low_rest_sex)      = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_sex))
colnames(design_tmt_low_rest_steroids) = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_steroids))
colnames(design_tmt_low_rest_therapy)  = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_therapy))
colnames(design_tmt_low_rest_epilepsy) = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_epilepsy))
colnames(design_tmt_low_rest_hg2a)     = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_hg2a))
colnames(design_tmt_low_rest_sept3)    = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_sept3))
colnames(design_tmt_low_rest_ptprc)    = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_ptprc))

# Base Model
colnames(design_tmt_low_rest_base)     = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_base))

# Full Model
colnames(design_tmt_low_rest_full)     = gsub("cond_tmt_low_rest", "", colnames(design_tmt_low_rest_full))





## Analysis  
# Limma 
fit_tmt_lr_1           <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest)

fit_tmt_lr_1_age       <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_age)
fit_tmt_lr_1_sex       <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_sex)
fit_tmt_lr_1_steroids  <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_steroids)
fit_tmt_lr_1_therapy   <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_therapy)
fit_tmt_lr_1_epilepsy  <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_epilepsy)
fit_tmt_lr_1_hg2a      <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_hg2a)
fit_tmt_lr_1_sept3     <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_sept3)
fit_tmt_lr_1_ptprc     <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_ptprc)

fit_tmt_lr_1_base      <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_base)

fit_tmt_lr_1_full      <- lmFit(dta_tmt_low_rest_conf[,-c(1:3)], design_tmt_low_rest_full)





# Contrast-Matrix of all comparisons no confounders
contrast_matrix_TMT <- makeContrasts(IDHHG-LowGrade, 
                                     RTKPN-LowGrade,
                                     RTKCL-LowGrade,
                                     MES-LowGrade,
                                     levels = design_tmt_low_rest)



# Univariable confouders
contrast_matrix_TMT_age      <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              age,
                                              levels = design_tmt_low_rest_age)

contrast_matrix_TMT_sex      <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              sexM,
                                              levels = design_tmt_low_rest_sex)

contrast_matrix_TMT_steroids <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              steroidsyes,
                                              levels = design_tmt_low_rest_steroids)

contrast_matrix_TMT_therapy  <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              therapy_anyTRUE,
                                              levels = design_tmt_low_rest_therapy)

contrast_matrix_TMT_epilepsy <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              epilepsyyes,
                                              levels = design_tmt_low_rest_epilepsy)

contrast_matrix_TMT_hg2a     <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              hg2a,
                                              levels = design_tmt_low_rest_hg2a)

contrast_matrix_TMT_sept3     <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              sept3,
                                              levels = design_tmt_low_rest_sept3)

contrast_matrix_TMT_ptprc     <- makeContrasts(IDHHG-LowGrade, 
                                              RTKPN-LowGrade,
                                              RTKCL-LowGrade,
                                              MES-LowGrade,
                                              ptprc,
                                              levels = design_tmt_low_rest_ptprc)


# Base
contrast_matrix_TMT_base <- makeContrasts(IDHHG-LowGrade, 
                                          RTKPN-LowGrade,
                                          RTKCL-LowGrade,
                                          MES-LowGrade,
                                          age,
                                          sexM,
                                          steroidsyes,
                                          therapy_anyTRUE,
                                          epilepsyyes,
                                          levels = design_tmt_low_rest_base)


# Full 

contrast_matrix_TMT_full <- makeContrasts(IDHHG-LowGrade, 
                                          RTKPN-LowGrade,
                                          RTKCL-LowGrade,
                                          MES-LowGrade,
                                          age,
                                          sexM,
                                          steroidsyes,
                                          therapy_anyTRUE,
                                          epilepsyyes,
                                          hg2a, 
                                          sept3,
                                          ptprc,
                                          levels = design_tmt_low_rest_full)



# Compute estimated coefficients and standard errors for a given set of contrasts
fit_tmt_lr_2      <- contrasts.fit(fit_tmt_lr_1, contrast_matrix_TMT)

# Univaribale counders
fit_tmt_lr_2_age      <- contrasts.fit(fit_tmt_lr_1_age, contrast_matrix_TMT_age)
fit_tmt_lr_2_sex      <- contrasts.fit(fit_tmt_lr_1_sex, contrast_matrix_TMT_sex)
fit_tmt_lr_2_steroids <- contrasts.fit(fit_tmt_lr_1_steroids, contrast_matrix_TMT_steroids)
fit_tmt_lr_2_therapy  <- contrasts.fit(fit_tmt_lr_1_therapy, contrast_matrix_TMT_therapy)
fit_tmt_lr_2_epilepsy <- contrasts.fit(fit_tmt_lr_1_epilepsy, contrast_matrix_TMT_epilepsy)
fit_tmt_lr_2_hg2a     <- contrasts.fit(fit_tmt_lr_1_hg2a, contrast_matrix_TMT_hg2a)
fit_tmt_lr_2_sept3    <- contrasts.fit(fit_tmt_lr_1_sept3, contrast_matrix_TMT_sept3)
fit_tmt_lr_2_ptprc     <- contrasts.fit(fit_tmt_lr_1_ptprc, contrast_matrix_TMT_ptprc)

fit_tmt_lr_2_base <- contrasts.fit(fit_tmt_lr_1_base, contrast_matrix_TMT_base)

fit_tmt_lr_2_full <- contrasts.fit(fit_tmt_lr_1_full, contrast_matrix_TMT_full)



# eBayes: Is estimating an "average" variability over all genes and is adjusting the individual variability accordingly 
# Uses all proteins
fit_tmt_lr_2      <- eBayes(fit_tmt_lr_2, trend = TRUE)

fit_tmt_lr_2_age       <- eBayes(fit_tmt_lr_2_age, trend = TRUE)
fit_tmt_lr_2_sex       <- eBayes(fit_tmt_lr_2_sex, trend = TRUE)
fit_tmt_lr_2_steroids  <- eBayes(fit_tmt_lr_2_steroids, trend = TRUE)
fit_tmt_lr_2_therapy   <- eBayes(fit_tmt_lr_2_therapy, trend = TRUE)
fit_tmt_lr_2_epilepsy  <- eBayes(fit_tmt_lr_2_epilepsy, trend = TRUE)
fit_tmt_lr_2_hg2a      <- eBayes(fit_tmt_lr_2_hg2a, trend = TRUE)
fit_tmt_lr_2_sept3     <- eBayes(fit_tmt_lr_2_sept3, trend = TRUE)
fit_tmt_lr_2_ptprc     <- eBayes(fit_tmt_lr_2_ptprc, trend = TRUE)

fit_tmt_lr_2_base <- eBayes(fit_tmt_lr_2_base, trend = TRUE)

fit_tmt_lr_2_full <- eBayes(fit_tmt_lr_2_full, trend = TRUE)



## DEPs
dep_IDHHG  <- topTable(fit_tmt_lr_2, coef = "IDHHG - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_RTKPN  <- topTable(fit_tmt_lr_2, coef = "RTKPN - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_RTKCL  <- topTable(fit_tmt_lr_2, coef = "RTKCL - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_MES    <- topTable(fit_tmt_lr_2, coef = "MES - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)


## Univariable Model 
dep_uni_age      <- topTable(fit_tmt_lr_2_age, coef = "age", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_uni_sex      <- topTable(fit_tmt_lr_2_sex, coef = "sexM", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_uni_steroids <- topTable(fit_tmt_lr_2_steroids, coef = "steroidsyes", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_uni_therapy  <- topTable(fit_tmt_lr_2_therapy, coef = "therapy_anyTRUE", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_uni_epilepsy <- topTable(fit_tmt_lr_2_epilepsy, coef = "epilepsyyes", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_uni_hg2a     <- topTable(fit_tmt_lr_2_hg2a, coef = "hg2a", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_uni_sept3    <- topTable(fit_tmt_lr_2_sept3, coef = "sept3", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_uni_ptprc    <- topTable(fit_tmt_lr_2_ptprc, coef = "ptprc", adjust.method="fdr", number=Inf, p.value = 0.05)    


##  Output full model
dep_base_IDHHG <- topTable(fit_tmt_lr_2_base, coef = "IDHHG - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_base_RTKPN <- topTable(fit_tmt_lr_2_base, coef = "RTKPN - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_base_RTKCL <- topTable(fit_tmt_lr_2_base, coef = "RTKCL - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_base_MES   <- topTable(fit_tmt_lr_2_base, coef = "MES - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)

dep_base_age      <- topTable(fit_tmt_lr_2_base, coef = "age", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_base_sex      <- topTable(fit_tmt_lr_2_base, coef = "sexM", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_base_steroids <- topTable(fit_tmt_lr_2_base, coef = "steroidsyes", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_base_therapy  <- topTable(fit_tmt_lr_2_base, coef = "therapy_anyTRUE", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_base_epilepsy <- topTable(fit_tmt_lr_2_base, coef = "epilepsyyes", adjust.method="fdr", number=Inf, p.value = 0.05)


##  Output full model
dep_full_IDHHG <- topTable(fit_tmt_lr_2_full, coef = "IDHHG - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_RTKPN <- topTable(fit_tmt_lr_2_full, coef = "RTKPN - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_RTKCL <- topTable(fit_tmt_lr_2_full, coef = "RTKCL - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_MES   <- topTable(fit_tmt_lr_2_full, coef = "MES - LowGrade", adjust.method="fdr", number=Inf, p.value = 0.05)

dep_full_age      <- topTable(fit_tmt_lr_2_full, coef = "age", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_sex      <- topTable(fit_tmt_lr_2_full, coef = "sexM", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_steroids <- topTable(fit_tmt_lr_2_full, coef = "steroidsyes", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_full_therapy  <- topTable(fit_tmt_lr_2_full, coef = "therapy_anyTRUE", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_full_epilepsy <- topTable(fit_tmt_lr_2_full, coef = "epilepsyyes", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_hg2a     <- topTable(fit_tmt_lr_2_full, coef = "hg2a", adjust.method="fdr", number=Inf, p.value = 0.05)
dep_full_sept3    <- topTable(fit_tmt_lr_2_full, coef = "sept3", adjust.method="fdr", number=Inf, p.value = 0.05)    
dep_full_ptprc    <- topTable(fit_tmt_lr_2_full, coef = "ptprc", adjust.method="fdr", number=Inf, p.value = 0.05)    




# --------------------------------
# Create Dataframe with all the DEP by contrast
# --------------------------------
for(i in list_of_comparison_low_rest){
    
    varname <- paste0("proteins_",gsub(" - ","_",i))
    assign(varname ,topTable(fit_tmt_lr_2_full, coef =  i, n = Inf,  p.value = 0.05, adjust.method = adjust_method))    
}

# List of DEP
tmt_lr_of_interest <- topTable(fit_tmt_lr_2_full, 
                               coef = c("IDHHG - LowGrade", "RTKPN - LowGrade", "RTKCL - LowGrade", "MES - LowGrade"),
                               n = Inf,  p.value = 0.05, 
                               adjust.method = adjust_method)



### Calculate the total number of significant DEP
## If we only consider p-values we get the same numbers as from the individual combination and number of p-values counted
results_lr = decideTests(fit_tmt_lr_2_full$p.value ,
                         method = "separate" ,
                         adjust.method = adjust_method)


n_decidetest_lr <- apply(results_lr[,1:4], 1, function(r) any(r == 1))
n_decidetest_lr <- n_decidetest_lr[n_decidetest_lr == TRUE]

## Names of all sinificant DEP overall contrasts
tmt_lr_of_interest <- attr(n_decidetest_lr,"names")

## List of DEP 
list_of_tmt_lr   <- paste0("proteins_",gsub(" - ","_", list_of_comparison_low_rest))
contrasts_TMT_lr <- gsub("proteins_", "" , list_of_tmt_lr)




# --------------------------------
# UP and DOWN regulated proteins
# --------------------------------

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
