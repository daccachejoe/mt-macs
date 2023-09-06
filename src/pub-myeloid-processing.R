# myeloid processing from the full annotated objects generated in pub-processing.R
myeloid.list <- list()
split.list <- list()
for(study in studies){
  load(file = paste0("data/new_analysis/", study, "_full_annotated.Rds"))
  myeloid.list[[study]] <- subset(object,
                                  cluster_annotation == "Myeloid" | cluster_annotation == "Cycling")
  rm(object)
  gc()
}

split.list <- lapply(myeloid.list, function(obj){
  new_list <- SplitObject(obj, split.by = "patient")
  return(new_list)
})

# Defining our integration function
# Find the conserved variable features as if doing Seurat CCA preprocess
HarmonyIntegration <- function(merged.object, split.list, vars.to.regress, FileName){
    
  split.list <- lapply(split.list, function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
  })
  
  int.features <- SelectIntegrationFeatures(object.list = split.list)
  message("Using n=", length(int.features), " to integrate.")
  
  # preprocess the objects using Seurat
  merged.object <- merged.object %>%
    NormalizeData(verbose = 3) %>%
    ScaleData(verbose = 3, features = int.features, vars.to.regress = vars.to.regress) %>%
    RunPCA(features = int.features, npcs = 30, verbose = 3)

  # RunHarmony
  merged.object <- merged.object %>%
    harmony::RunHarmony("patient", 
                        plot_convergence = FALSE, 
                        max.iter.harmony = 50)

  merged.object <-  FindNeighbors(merged.object, 
                                  reduction = "harmony", 
                                  graph.name = c("harmony_nn","harmony_snn"),
                                  dims = 1:30)

    # we ran this to know what resolutions to use, then commented it out
  # resolutions_to_test <- seq(0,1,.1)
  # for(res in resolutions_to_test){
  #   merged.object <- FindClusters(merged.object, graph.name = "harmony_snn", resolution = resolutions_to_test)
  # }
  
  red_to_use <- "umap"
  res_to_use <- list(lung = 0.5,
                     liver = 0.6,
                     kidney = 0.4)
  merged.object <- FindClusters(merged.object, graph.name = "harmony_snn", resolution = res_to_use[[merged.object$organ[1]]])
  merged.object <- RunUMAP(merged.object,dims = 1:30, reduction = "harmony")
  return(merged.object)
}

# integrating the Myeloid cells in each study
for(study in studies){
  gc()
  myeloid.list[[study]] <- 
    HarmonyIntegration(merged.object = myeloid.list[[study]],
                       split.list = split.list[[study]],
                       vars.to.regress = "patient",
                       FileName = paste0("data/new_analysis/", 
                                         study,
                                         "_myeloid.Rds"))
  gc()
}

names(myeloid.list) <- c("lung", "liver", "kidney")
myeloid.list <- lapply(myeloid.list, function(obj){
  obj@project.name <- obj$organ[1]
  return(obj)
})


# selecting resolutions to continue with
library(clustree)
clus.tree.plots <- lapply(myeloid.list, function(mye){
  p <- clustree::clustree(mye,layout = "tree", prefix = "harmony_snn_res.")
  return(p)
})
clus.tree.plots

dim.plot.list <- lapply(myeloid.list, function(mye){
  factors_to_plot <- paste0("harmony_snn_res.", seq(0,1,.1))
  p.list <- lapply(factors_to_plot, function(fac){
    p <- 
      DimPlot(mye, group.by = fac,
              pt.size = 0.1, label = T) + 
      NoLegend() +
      theme(axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()) +
      ggtitle(fac)
    return(p)
  })
  return(p.list)
})
ggpubr::ggarrange(plotlist = dim.plot.list[[1]], nrow = 3, ncol = 4)
ggpubr::ggarrange(plotlist = dim.plot.list[[2]], nrow = 3, ncol = 4)
ggpubr::ggarrange(plotlist = dim.plot.list[[3]], nrow = 3, ncol = 4)

# start of interactive and iterative analysis 
# now we have selected a resolution and clustering result to use
red_to_use <- "umap"
res_to_use <- list(lung = "harmony_snn_res.0.5",
                   liver = "harmony_snn_res.0.6",
                   kidney = "harmony_snn_res.0.4")

# UMAP plots of new cluster definitions
dimplot.list <- lapply(myeloid.list, function(obj){
  p <- DimPlot(obj, group.by = res_to_use[[obj@project.name]],
               label = T, pt.size = 0.01) + NoLegend() + ggtitle(obj@project.name)
  return(p)
})
vlnplot.list <- lapply(myeloid.list, function(obj){
  p <- VlnPlot(obj, group.by = res_to_use[[obj@project.name]],
               pt.size = 0, 
               features = c("S100A8", "FCN1", "C1QB","MARCO", "CD163","LYVE1",
                            "TREM2", "CD9", "SPP1", "RNASE1","FOLR2", "CXCL10","MT1G",
                            ifelse(obj@project.name == ("liver", "SEPP1", "SELENOP"),
                            "LIPA", "LPL", "CD36", "FABP4", "FABP5"), stack = T) +
    ggtitle(obj@project.name) +
    NoLegend()
})
ggpubr::ggarrange(plotlist = dimplot.list, nrow = 2, ncol = 2)
vlnplot.list

# Find the cluster markers of the clusters
marker.list <- lapply(myeloid.list, function(obj){
  Idents(obj) <- res_to_use[[obj@project.name]]
  M <- FindAllMarkers(obj, only.pos = T, max.cells.per.ident = 5000)
  return(M)
})
writexl::write_xlsx(x = marker.list, path = "output/myeloid_cluster_markers.xlsx")


# read in our annotations, remove problematic clusters, then re-run clustering with only monocytes and macrophages
annotations <- lapply(c("lung", "liver", "kidney"), function(x){
  M <- read.csv(paste0("output/", x, "_myeloid_cluster_idens.csv"), col.names = c("cluster", "myeloid_cellType", "group_var"))
  return(M)
})
names(annotations) <- c("lung", "liver", "kidney")

# major macrophage levels defined here
mac_levels <- read.csv("output/releveling_macrophages.csv")

for(n in names(myeloid.list)){
  md.to.add <- data.frame(row.names = rownames(myeloid.list[[n]]@meta.data),
                          myeloid_cellType = annotations[[n]]$myeloid_cellType[match(myeloid.list[[n]]@meta.data[ , res_to_use[[n]]], 
                                                                                     annotations[[n]]$cluster)],
                          group_var = annotations[[n]]$group_var[match(myeloid.list[[n]]@meta.data[ , res_to_use[[n]]],
                                                                                      annotations[[n]]$cluster)])
  
  myeloid.list[[n]] <- AddMetaData(myeloid.list[[n]], md.to.add)
}

# remove identified provlematic clusters
for(n in names(myeloid.list)){
  myeloid.list[[n]] <- subset(myeloid.list[[n]], group_var %in% mac_levels$group_var)
}

# save the objects
for(object in myeloid.list){
  save(object, file = paste0("data/new_analysis/", object@project.name, "_annotated_mono_mac_prolif.Rds"))
}

# Combine everything into one object and move onto figure panels script
object.list <- lapply(c(
    "liver",
    "kidney",
    "lung"
), function(organ) {
    author <- ifelse(organ == "lung", "Adams_Ayub",
        ifelse(organ == "liver",
            "Ramachandran",
            "Kuppe"
        )
    )
    full.object.path <- paste0("r-objects/", author, "_full_annotated.Rds")
    load(file = full.object.path)
    myeloid.path <- paste0("r-objects/", organ, "_annotated_mono_mac_prolif.Rds")
    load(file = myeloid.path)
    myeloid.object <- object

    md.to.add <- data.frame(
        row.names = rownames(myeloid.object@meta.data),
        group_var = myeloid.object$group_var
    )
    md.to.add$group_var <- factor(md.to.add$group_var, levels = unique(mac_levels$group_var))
    myeloid.object$group_var <- md.to.add$group_var

    rm(object)
    gc()

    # bringing over the myeloid annotations to the full object
    indepth_CT_table <-
        data.frame(
            cells = colnames(full.object),
            indepth_cellType = full.object$cluster_annotation
        )
    indepth_CT_table$indepth_cellType[
        which(indepth_CT_table$cells %in% colnames(myeloid.object))
    ] <-
        as.character(myeloid.object$group_var)
    full.object <-
        AddMetaData(full.object, metadata = indepth_CT_table)

        # remove previously identified myeloid lineage cells that didnt pass the sub-clustering filtering part
    full.object$ct.filtered <- ifelse(full.object$indepth_cellType == "Myeloid" &
        !(full.object$individual_CT_Specific %in% c("cDC 1", "cDC 2")),
    "rm", "keep"
    )
    full.object <- subset(full.object, ct.filtered == "keep")
    full.object$indepth_cellType <- ifelse(full.object$indepth_cellType == "Myeloid",
        full.object$individual_CT_Specific,
        full.object$indepth_cellType
    )
   
    full.object <- RunUMAP(full.object,
        reduction = "harmony",
        dims = 1:50
    )
    return(list(
        full = full.object,
        myeloid = myeloid.object
    ))
})
names(object.list) <- c(
    "liver",
    "kidney",
    "lung"
)
save(object.list, file = "r-objects/object.list.rds")