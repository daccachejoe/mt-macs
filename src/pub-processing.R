# dimensionality reduction pipeline run on each object then subsetted down to macrophages only


ProcessingPipeline <- function(object, projectName, Cores){
    if(projectName == "Adams/Ayaub"){
        projectName <- "lung"
    }
    IntegratedObject <- IntegrateMyData(object,
                            projectName = projectName,
                            dims = 50,
                            assay = "RNA",
                            method = "harmony",
                            HarmonyVar =c("sample"),
                            harmonyTheta = c(2),
                            cores = Cores)
    save(IntegratedObject, file = paste0("data/combined/Seurat/", projectName, ".Rds"))

    pdf(file = paste0("output/qc/PCA_before_after_", projectName, ".pdf"), height = 8, width = 10)
        p1 <- DimPlot(IntegratedObject, reduction = "pca", pt.size = 0.01, group.by = "organ") + ggtitle("pre-correction") + BoldTitle()
        p2 <- DimPlot(IntegratedObject, reduction = "harmony", pt.size = 0.01, group.by = "organ") + ggtitle("post-correction") + BoldTitle()
        print(ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "top"))
    dev.off()

    # Running the determined Dimensionality Reduction pipeline with heightened dimensions, resolutions, and perplexities given large cell counts
    message("Running Dimensionality Reduction")
    IntegratedObject <- RunDimensionalityReduction(object = list(IntegratedObject),
                                                    assay = "RNA",
                                                    dims = 50, 
                                                    cores = Cores,
                                                    vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt","S.Score", "G2M.Score"),
                                                    resolution = c(.5, 1),
                                                    perplexity = c(50, 100),
                                                    use.reduction = "harmony")
    message("Dimensionality Reduction complete")
    save(IntegratedObject, file = paste0("data/combined/Seurat/", projectName, ".Rds"))

    message("Processing pipeline complete on ", projectName)
    return(IntegratedObject)
}

# Bring in only protein coding genes
# HS_gene_info <- annot[annot$`#tax_id` == 9606, ]
# write.table(HS_gene_info, file = "data/homo_sapien_gene_info.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# Start by consolidating gene names
annot <- readr::read_delim("data/homo_sapien_gene_info.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
conversion_df <-
  annot %>% 
  filter(type_of_gene == "protein-coding")
PC_Genes <- conversion_df$Symbol_from_nomenclature_authority
Syn_Genes <- unlist(strsplit(conversion_df$Synonyms, split = "|", fixed = TRUE))
Full_Genes <- c(PC_Genes, Syn_Genes)

# export clusters and author annotations
for(study in studies){
    # Load in the R objects
    load(file = paste0("data/", study, "_Annotated.Rds"))
    gc()
  
    # Exporting the author lineages percentage by RNA cluster and our previous annotations (not shown, iteratively coded)
    df <- as.data.frame(table(object$RNA_snn_res.1, 
                              object$cellType_Specific))
    colnames(df) <- c("cluster", "author_annot", "Freq")
    df <- df %>%
      filter(Freq != 0) %>%
      group_by(author_annot) %>%
      mutate(perc = Freq/sum(Freq))
    write.csv(df, file = paste0("output/", study, "_ct_table.csv")) 
    
    # Transferring author lineage labels to our clusters
    df <- as.data.frame(table(object$RNA_snn_res.1, 
                              object$cell_Lineage))
    colnames(df) <- c("cluster", "author_annot", "Freq")
    df <- df %>%
      filter(Freq != 0) %>%
      group_by(cluster) %>%
      mutate(perc = Freq/sum(Freq))
    new_annotations <- data.frame(cluster = unique(df$cluster),
                                  new_annotation = character(length = length(unique(df$cluster))))
    for(my_cluster in unique(df$cluster)){
      df_to_retain <- df[df$cluster == my_cluster, ]
      if(max(df_to_retain$perc) > 0.8){
        label <- df_to_retain$author_annot[which.max(df_to_retain$perc)]
        message("Cluster ", my_cluster, " annotated as ", label)
        new_annotations$new_annotation[new_annotations$cluster == my_cluster] <- as.character(label)
      } else {
        label <- " "
        message("Cluster ", my_cluster, " requires manual annotation.")
        new_annotations$new_annotation[new_annotations$cluster == my_cluster] <- label
      }
    }
    write.csv(new_annotations,
              row.names = FALSE,
              file = paste0("output/", study, "_author_lineages_converted_to_daccache_clusters.csv"))
    
    # clean-up
    rm(object)
    gc()
  }

  # After manual inputs are made as needed, bring the annotations back into R and add them to the objects
  if (ReCreateSeuratObject) {
      # After manual annotations are complete, read back in the data frames
      annotations <- lapply(studies, function(study) {
          df <- read.csv(file = paste0("output/", study, "_author_lineages_converted_to_daccache_clusters.csv"))
          return(df)
      })
      names(annotations) <- studies

      # Cross-study annotations
      conserved_lineage_names <- read.csv(file = "output/conserving_lineage_labels.csv", col.names = c("converted_lineage_name", "author_annotation"))

      # Now we add them one by one
      for (study in studies) {
          message("Working on ", study)
          load(file = paste0("data/", study, "_Annotated.Rds"))

          gc()

          # Quick pivot for the lung study
          if (study == "Adams_Ayub") {
              object$author <- "Adams_Ayub"
          }

          # Re-formatting the meta data DF's to contain the columns we want to keep and in a legible order
          cols_to_keep <- c(
              "orig.ident", "hash", "barcode", "BC_ranked",
              "nCount_RNA", "nFeature_RNA",
              "study", "author", "organ", "tissue", "patient",
              "sample", "sex", "age", "Pre_Ex", "collection",
              "condition", "cause_of_disease",
              "S.Score", "G2M.Score", "Phase",
              "cell_Lineage", "cellType", "cellType_Specific",
              "RNA_snn_res.0.5", "RNA_snn_res.1", "seurat_clusters"
          )
          md <- object@meta.data[, cols_to_keep]
          message("adding annotations")

          annotations[[unique(object$author)[1]]]$cluster <- as.character(annotations[[unique(object$author)[1]]]$cluster)
          md <- left_join(md, annotations[[unique(object$author)[1]]], by = c("RNA_snn_res.1" = "cluster"))
          colnames(md)[which(colnames(md) == "new_annotation")] <- "cluster_annotation"

          rownames(md) <- md$barcode
          object@meta.data <- md

          # This is the tough one:
          # Removal of dual-lineage barcodes and non PC genes
          # NOTE: have to run lung study by itself unless have a CPU with memory > 50G
          # Did not run with large memory environment on desktop with 50G of RAM
          # not sure what memory will work
          message("Cross-referencing to published annotations")
          reference_df <- conserved_lineage_names[conserved_lineage_names$author_annotation %in% unique(md$cell_Lineage), ]
          reference_df <- rbind(reference_df, conserved_lineage_names[conserved_lineage_names$author_annotation %in% c("keep", "remove"), ])

          md$converted_author_annotation <- reference_df$converted_lineage_name[match(md$cell_Lineage, reference_df$author_annotation)]
          md$to_keep <- if_else(md$converted_author_annotation == "keep", TRUE,
              ifelse(md$converted_author_annotation == "remove", FALSE,
                  ifelse(md$cluster_annotation == "Plasma", TRUE,
                      ifelse(md$cluster_annotation == "Cycling", TRUE,
                          ifelse(md$converted_author_annotation == md$cluster_annotation,
                              TRUE,
                              FALSE
                          )
                      )
                  )
              )
          )
          if (study == "Kuppe") {
              md$to_keep <- ifelse(md$to_keep == FALSE,
                  ifelse(md$cell_Lineage == "Immune", TRUE, FALSE),
                  md$to_keep
              )
          }
          remove_list <- as.data.frame(table(md$cluster_annotation, md$cell_Lineage, md$to_keep))
          colnames(remove_list) <- c("cluster_annotation", "author_annotation", "keep", "counts")
          write.csv(remove_list,
              file = paste0(
                  "output/",
                  unique(object$author),
                  "_converted_lineage_to_remove.csv"
              ),
              row.names = FALSE
          )

          errored_groups <-
              remove_list %>%
              group_by(author_annotation, keep) %>%
              summarise(kept_counts = sum(counts)) %>%
              mutate(total_counts = sum(kept_counts)) %>%
              mutate(flag = ifelse(kept_counts == total_counts, TRUE, FALSE)) %>%
              filter(flag) %>%
              filter(keep != "TRUE")

          if (nrow(errored_groups) > 0) {
              if (unique(errored_groups$author_annotation) != "Neuronal" & unique(errored_groups$author_annotation) != "Multiplet") { # we want all Neuronal kidney and Multiplet lung cells removed
                  stop("author groups: ", paste(unique(errored_groups$author_annotation), collapse = " "), " all removed")
              }
          }

          
          if (study == "Adams_Ayub") {
              md$to_keep <- ifelse(md$cause_of_disease == "COPD", FALSE, md$to_keep)
          }

          cells_to_keep <- rownames(md)[md$to_keep]
          md <- md[cells_to_keep, ]

          message("Trimming down Seurat object")

          # Consolidating gene names via synonyms
          # Retaining only PC genes
          old.matrix <- object[["RNA"]]@counts[rownames(object[["RNA"]]@counts) %in% Full_Genes, ]
          dim_reducs_to_keep <- object@reductions[c("pca", "harmony", "RNA_tsne_100")]
          rm(object)
          gc()
          gc()

          # filtering of matrix for genes expressed in at least 10 cells
          genes_to_keep <- rownames(old.matrix)[rowSums(old.matrix > 0) > 10]
          # creating metric that takes total counts/total # of positive cells
          gene_depth <- rowSums(old.matrix[genes_to_keep, ]) / rowSums(old.matrix[genes_to_keep, ] > 0)
          genes_to_keep <- names(gene_depth)[gene_depth > 1]
          # remove features with a 1:1 ratio of positive cells and total roads (i.e. 1 read/cell)
          old.matrix <- old.matrix[genes_to_keep, ]
          gc()

          object <- CreateSeuratObject(
              counts = old.matrix,
              meta.data = md
          )
          object@project.name <- study
          object@reductions <- dim_reducs_to_keep

          message("Processing complete on ", study)
          save(object, file = paste0("data/new_analysis/", study, "_full_annotated.Rds"))
          gc()
      }
  }
