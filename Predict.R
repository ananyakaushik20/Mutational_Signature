#load dependencies
library(dplyr)
library(Matrix)
library(knitr)
library(irlba)

#load required packages
library(renv)
library(TCGAbiolinks)
library(withr)
library(musicatk)
library(BSgenome)

# Import TCGAbiolinks data----
tcga_datasets <- c("TCGA-COAD","TCGA-READ")

types <- sapply(strsplit(tcga_datasets, "-"), '[', 2)
dataset <- NULL
annot <- NULL
for (i in seq_along(tcga_datasets)){
  query <- GDCquery(project=tcga_datasets[i],
                    data.category = "Simple Nucleotide Variation",
                    data.type = "Masked Somatic Mutation",
                    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                    experimental.strategy = "WXS",
                    data.format = "maf")
  GDCdownload(query)
  data <- GDCprepare(query)
  variants <- extract_variants_from_matrix(data)
  dataset <- rbind(dataset, variants)
  annot <- rbind(annot, cbind(rep(types[i], length(unique(variants$sample))),
                              unique(as.character(variants$sample))))
}

#create musica object
g <- select_genome("hg38")
tcga <- create_musica(dataset,g)

#sample annotation
matched_annot <- annot[match(samp_annot(tcga)$Samples,annot[,2]),1]
samp_annot(tcga, "Tumor_Types") <- matched_annot

# SBS Umap

build_standard_table(tcga,g, table_name = "SBS96")

#filtering low count samples
tcga_subset <- subset_musica_by_counts(tcga, table_name = "SBS96", num_counts = 5)

tcga_v3_sbs <- with_seed(1, auto_predict_grid(tcga_subset, 
                                              table_name = "SBS96",
                                              signature_res = cosmic_v3_sbs_sigs_exome,
                                              algorithm = "lda",
                                              sample_annotation = "Tumor_Types")
  )

plot_exposures(tcga_v3_sbs, 
               group_by = "annotation", 
               annotation = "Tumor_Types",
               proportional = TRUE)

with_seed(1,create_umap(result = tcga_v3_sbs))

plot_umap(result = tcga_v3_sbs, 
          proportional = TRUE, 
          color_by = "annotation",
          annotation = "Tumor_Types",
          add_annotation_labels = TRUE,
          annotation_text_box = TRUE,
          annotation_label_size = 6,
          legend = FALSE,
          strip_axes = TRUE)

plot_umap(result = tcga_v3_sbs,same_scale = FALSE)
