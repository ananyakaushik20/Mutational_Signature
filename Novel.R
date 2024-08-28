#Novel Signatures


denovo <- discover_signatures(tcga_subset, "SBS96", num_signatures = 10, nstart = 1)

plot_signatures(denovo)

compare_cosmic_v2(denovo)

plot_exposures(denovo,
               proportional = TRUE, 
               sort_samples = "Signature9", 
               group_by = "annotation", 
               annotation = "Tumor_Types")

plot_heatmap(res_annot = denovo, proportional = TRUE, scale = TRUE, annotation = "Tumor_Types")

#Linear Regression

glm <- exposure_differential_analysis(denovo, 
                                      annotation = "Tumor_Types",
                                      method = "glm.nb")

plot_differential_analysis(glm, "glm", 2)