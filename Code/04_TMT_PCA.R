# Function to replace a pattern in a vector of strings
replace_string<- function(string_vector, pattern, replacement) {
    # Check if the inputs are valid
    if(!is.vector(string_vector) || !is.character(string_vector)) {
        stop("Input 'string_vector' must be a character vector.")
    }
    
    if(!is.character(pattern) || length(pattern) != 1) {
        stop("Input 'pattern' must be a single character string.")
    }
    
    if(!is.character(replacement) || length(replacement) != 1) {
        stop("Input 'replacement' must be a single character string.")
    }
    
    # Replace the pattern in each string of the vector
    replaced_vector <- gsub(pattern, replacement, string_vector)
    
    return(replaced_vector)
}



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
colnames(dta_tmt_PCA) <- replace_string(colnames(dta_tmt_PCA), "_HUMAN", "")



# --------------------------------------------
# PCA
# --------------------------------------------

res_pca   <- PCA(dta_tmt_PCA[, !(names(dta_tmt_PCA) %in% "groups")],  graph = FALSE)

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
# get_eig(res_pca)
# fviz_eig(res_pca)
fviz_screeplot(res_pca, addlabels = TRUE, ylim = c(0, 40))

# Graph of variables. Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(res_pca, 
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
             )


# Biplot of individuals and variables
fviz_pca_biplot(res_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )

var <- get_pca_var(res_pca)
corrplot(var$cos2, is.corr = FALSE)

# Contribution of each variable 
fviz_cos2(res_pca, 
          choice = "var", 
          axes = 1:2, # Which Dimensions of the PCA
          top = 50)   # Top 50 variables


# Graph of individuals. Individuals with a similar profile are grouped together.
fviz_pca_ind(res_pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


pca_DEP <- fviz_pca_ind(res_pca,
             palette = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
             label = "none", # hide individual labels
             habillage = dta_tmt_PCA$groups, # color by groups
             # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             title = "PCA visualization"
) + theme(legend.title=element_blank()) + theme(legend.justification = "top") +
    theme(
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18)  # center title
    )
# + guides(fill = "none")


pca_DEP_2 <- fviz_pca_biplot(res_pca, repel = FALSE,
                # col.var = "#2E9FDF", # Variables color
                col.var = "#999999", # Variables color
                alpha.var = 0.5,
                # alpha.ind = 0.1,
                # col.ind = "#696969",  # Individuals color
                palette = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
                label = "var", # hide individual labels
                habillage = dta_tmt_PCA$groups, # color by groups
                # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                addEllipses = TRUE, # Concentration ellipses
                title = "PCA visualization"
)


