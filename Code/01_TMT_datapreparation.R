###########################
## Project: Proteomics
## Data preparation 
## Created: 17.09.2020
###########################

## Example article
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8009917/


# Manual exclude patients
dta_full_tmt <- dta_full_tmt[,!(colnames(dta_full_tmt) %in% "IDHTERT_027_exp1")]
dta_full_tmt <- dta_full_tmt[,!(colnames(dta_full_tmt) %in% "RTKCL_061_exp7")]
dta_full_tmt <- dta_full_tmt[,!(colnames(dta_full_tmt) %in% "MES_072_exp7")]
dta_full_tmt <- dta_full_tmt[,!(colnames(dta_full_tmt) %in% "MES_073_exp7")]

### Data preparation
# Are there any duplicated gene names?
dta_full_tmt$UniprotID %>% duplicated() %>% any()
colnames(dta_full_tmt)[colnames(dta_full_tmt) == "UniprotID"]     <- "ID"
colnames(dta_full_tmt)[colnames(dta_full_tmt) == "Protein.Names"] <- "name"

## Cleaning protein names
rows_semicolcon <- grep(";", dta_full_tmt$ID)
dta_full_tmt$ID[rows_semicolcon]

## Cleaning-up protein names 
## Consider only first protein name
dta_full_tmt$ID[rows_semicolcon] <- str_split(dta_full_tmt$ID[rows_semicolcon], ";", simplify = T)[,1]


#######################################################
## TMT Proteomics
#######################################################
## Data preparation 
# Keep only abundance measures 
snp_cols <- grep("^SNP_", names(dta_full_tmt))
dta_tmt <- dta_full_tmt[,-snp_cols]

# dta_tmt <- dta_full_tmt[,-c(90:99)]
row.names(dta_tmt) <- dta_tmt$name

## Log transform abundance measures 
dta_tmt_log <- dta_tmt
dta_tmt_log[,4:85] <-  log2(dta_tmt_log[,4:85])
rownames(dta_tmt_log) <- dta_tmt_log$name

### DATA preparation
# Filtering for NAs
# Cut-off of proportion of allowed missings
# If cut-off = 1, only complete cases will kept.
cutoff_propgroup <- 1

# Missings by sample 
# Returns proportion of samples which are not missing for each protein and group
cols_IDHLG <- grep("IDHLG",colnames(dta_tmt_log))
keep_na_IDHLG   <- apply(dta_tmt_log[,cols_IDHLG],  1, 
                       function(x) sum(!is.na(x))/ length(cols_IDHLG) )

cols_IDHHG <- grep("IDHHG",colnames(dta_tmt_log))
keep_na_IDHHG   <- apply(dta_tmt_log[,cols_IDHHG],  1, 
                       function(x) sum(!is.na(x))/ length(cols_IDHHG))

cols_IDHTERT <- grep("IDHTERT",colnames(dta_tmt_log))
keep_na_IDHTERT   <- apply(dta_tmt_log[,cols_IDHTERT],  1, 
                      function(x) sum(!is.na(x))/ length(cols_IDHTERT))

cols_IDHCODEL <- grep("IDHCODEL",colnames(dta_tmt_log))
keep_na_IDHCODEL   <- apply(dta_tmt_log[,cols_IDHCODEL],  1, 
                       function(x) sum(!is.na(x))/ length(cols_IDHCODEL))

cols_RTKPN <- grep("RTKPN",colnames(dta_tmt_log))
keep_na_RTKPN   <- apply(dta_tmt_log[,cols_RTKPN],  1, 
                       function(x) sum(!is.na(x))/ length(cols_RTKPN))

cols_RTKCL <- grep("RTKCL",colnames(dta_tmt_log))
keep_na_RTKCL   <- apply(dta_tmt_log[,cols_RTKCL],  1, 
                       function(x) sum(!is.na(x))/ length(cols_RTKCL))

cols_MES <- grep("MES",colnames(dta_tmt_log))
keep_na_MES   <- apply(dta_tmt_log[,cols_MES],  1, 
                        function(x) sum(!is.na(x))/ length(cols_MES))

# Keep only the proteins where each group has more samples than the proposed cut-off
keep_tmt_na <- 
  ((keep_na_IDHLG   >= cutoff_propgroup) + 
  (keep_na_IDHHG    >= cutoff_propgroup) + 
  (keep_na_IDHTERT  >= cutoff_propgroup) + 
  (keep_na_IDHCODEL >= cutoff_propgroup) + 
  (keep_na_RTKPN    >= cutoff_propgroup) + 
  (keep_na_RTKCL    >= cutoff_propgroup) + 
  (keep_na_MES      >= cutoff_propgroup))  >= 7

# Exclude Proteins with more missing samples than the defined cutoff number
dta_tmt_log <- dta_tmt_log[keep_tmt_na,]





# ------------------------------
# Counfounder data preparation
# ------------------------------

# Name of colums - for group names 
names_cols <- sub("(_exp.*)", "", colnames(dta_tmt_low_rest_conf[,-c(1:3)]))

# Rename groups to fit to TMT data set
dta_confounder$Subject <- str_replace(dta_confounder$Subject, "GBMES", "MES")
dta_confounder$Subject <- str_replace(dta_confounder$Subject, "GBPN", "RTKPN")
dta_confounder$Subject <- str_replace(dta_confounder$Subject, "GBCL", "RTKCL")
dta_confounder$Subject <- str_replace(dta_confounder$Subject, "IDHLGG", "IDHLG")
dta_confounder$Subject <- str_replace(dta_confounder$Subject, "IDHHGG", "IDHHG")

# Patients with id more than 25 were called IDHTERT in the TMT data set
nums <- as.numeric(sub("IDHLG_", "", dta_confounder$Subject[1:23]))
dta_confounder$Subject[1:23] <- ifelse(nums >= 25, sub("IDHLG", "IDHTERT", dta_confounder$Subject[1:23]), dta_confounder$Subject[1:23])

dta_confounder$therapy_any <- dta_confounder$Chemotherapy == "yes" | dta_confounder$Radiotherapy == "yes"
dta_confounder <- dta_confounder[match(names_cols, dta_confounder$Subject), ]

# Load confounders from TMT data set
var_sept3    <- dta_tmt_low_rest_conf[dta_tmt_low_rest_conf$name == "SEPT3_HUMAN",]
var_sept3    <- as.numeric(var_sept3[,-c(1:3)])
var_ptprc    <- dta_tmt_low_rest_conf[dta_tmt_low_rest_conf$name == "PTPRC_HUMAN",]
var_ptprc    <- as.numeric(var_ptprc[,-c(1:3)])


age         <- as.numeric(dta_confounder$`Age OP`)
sex         <- as.factor(dta_confounder$Gender)
steroids    <- as.factor(dta_confounder$Steroids)
therapy_any <- as.factor(dta_confounder$therapy_any)
epilepsy    <- as.factor(dta_confounder$Epilepsy)
hg2a        <- log2(dta_confounder$HG2A_HUMAN)
sept3       <- var_sept3
ptprc       <- var_ptprc
