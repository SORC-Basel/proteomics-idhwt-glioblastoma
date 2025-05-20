

## Vulcano plot
# This function is saved in the separated custom function files
p_vulcano_idhhg <- ggplot_vulcano(fit_tmt_lr_2_full, "IDHHG - LowGrade", p_value = "unadjusted", labels = TRUE )
p_vulcano_rtkpn <- ggplot_vulcano(fit_tmt_lr_2_full, "RTKPN - LowGrade", p_value = "unadjusted", labels = TRUE )
p_vulcano_rtkcl <- ggplot_vulcano(fit_tmt_lr_2_full, "RTKCL - LowGrade", p_value = "unadjusted", labels = TRUE )
p_vulcano_mes <- ggplot_vulcano(fit_tmt_lr_2_full, "MES - LowGrade", p_value = "unadjusted", labels = TRUE )




