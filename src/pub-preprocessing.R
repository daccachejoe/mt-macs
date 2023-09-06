# preprocessing of data from published sources
# harmonizing cell level meta data 

Cores <- as.numeric(args[1]) # please pass the number of alloted cores in the bash script after the R script is called
library(future)
library(doFuture)
library(Seurat)
library(dplyr)
library(Matrix)
registerDoFuture()
plan("multiprocess", workers = Cores)

# (1) Preprocessing the lung IPF/COPD/Healthy data
# read into the procesed matrices from the lung data
message("Reading in lung data")
Matrix <- read.table("data/lung/outs/Pub_Raw_counts.txt", header = FALSE, skip = 2)
Features <- read.table("data/lung/outs/features.tsv.gz", header = FALSE)
Barcodes <- read.table("data/lung/outs/barcodes.tsv", header = FALSE)
counts <- Matrix::sparseMatrix(i = Matrix$V1, j = Matrix$V2, x = Matrix$V3)
rownames(counts) <- Features$V2
colnames(counts) <- Barcodes$V1

# read in the meta data
message("Adding meta data to merged lung object")
Pub_Meta_Data <- read.table("data/lung/Seurat/GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", header = TRUE)
Patient_MetaData <- read.csv("output/qc/lung_patient_ID_and_MetaData.csv", header = TRUE)
Pub_Meta_Data <- subset(Pub_Meta_Data, CellBarcode_Identity %in% colnames(counts))
Pub_Meta_Data <- left_join(Pub_Meta_Data, Patient_MetaData, by = c("Subject_Identity" = "PID"))
row.names(Pub_Meta_Data) <- Pub_Meta_Data$CellBarcode_Identity
Merged_Lung <- CreateSeuratObject(counts, min.cells = 3, meta.data = Pub_Meta_Data)
Merged_Lung <- CellCycleScoring(Merged_Lung, s.features = toupper(cc.genes$s.genes), toupper(cc.genes$g2m.genes))

# rename meta data columns to harmonize the column names
message("Re-naming meta data columns")
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "CellBarcode_Identity")] <- "hash"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "CellType_Category")] <- "cell_Lineage"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "Manuscript_Identity")] <- "cellType"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "Subclass_Cell_Identity")] <- "cellType_Specific"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "Disease_Identity")] <- "cause_of_disease"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "Subject_Identity")] <- "patient"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "Library_Identity")] <- "sample"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "nUMI")] <- "nCount_RNA"
colnames(Merged_Lung@meta.data)[which(colnames(Merged_Lung@meta.data) == "nGene")] <- "nFeature_RNA"
Merged_Lung$condition <- ifelse(Merged_Lung$cause_of_disease == "Control", "Uninjured", "Diseased")
Merged_Lung[["percent.mt"]] <- PercentageFeatureSet(Merged_Lung, pattern = "^MT-")
Merged_Lung[["percent.ribo"]] <- PercentageFeatureSet(Merged_Lung, pattern = "^RP")
Merged_Lung$orig.ident <- Merged_Lung$sample
Merged_Lung$organ <- "lung"
Merged_Lung$collection <- "sc_dissociation"
Merged_Lung$study <- "Adams_Ayaub"
Merged_Lung$tissue <- "lung_explant"
Merged_Lung$barcode <- sapply(strsplit(Merged_Lung$hash, split = "_"), getElement, 2)

message("Saving merged lung object")
save(Merged_Lung, file = "data/lung/Seurat/Merged_Lung.Rds")

# (2) Now we load in the liver data from Ramachandran 2019
message("Reading in liver data")
# read in the 10x matrices from the liver data
LiverFiles <- list.dirs("data/Ramachandran/10x_output", recursive=FALSE, full.names = FALSE)
Liver_List <- lapply(LiverFiles, function(x){
    filePath <- paste0("data/Ramachandran/10x_output/", x, "/outs/filtered_feature_bc_matrix")
    counts <- Read10X(filePath)
    SeuObject <- CreateSeuratObject(counts, min.cells = 0, min.features = 0, project = x)
})
message("Merging liver runs")
Merged_Liver <- merge(x = Liver_List[[1]], y = Liver_List[c(2:length(Liver_List))])
Merged_Liver <- CellCycleScoring(Merged_Liver, s.features = toupper(cc.genes$s.genes), toupper(cc.genes$g2m.genes))
Rama_samplesheet <- read.csv("../Human_Liver_sc_Comparison/samplesheet_Ramachandran.csv")
Merged_Liver@meta.data <- left_join(Merged_Liver@meta.data, Rama_samplesheet, by = c("orig.ident" = "projectName"))
Merged_Liver$barcode <- sapply(strsplit(colnames(Merged_Liver), split = "-"), getElement, 1)
Merged_Liver$hash <- paste(Merged_Liver$Sample, Merged_Liver$barcode, sep = "_")

message("Adding meta data to merged liver object")
# read in the published meta data from the liver data and add it to the 10x matrices we read in above
load("data/Ramachandran/published_Seurat/tissue.rdata")
Pub_Seurat <- CreateSeuratObject(counts = tissue@data, meta.data = tissue@meta.data, min.cells = 0, min.features = 0)
Pub_Seurat$hash <- colnames(Pub_Seurat)

rownames(Merged_Liver@meta.data) <- colnames(Merged_Liver)
Merged_Liver <- subset(Merged_Liver, hash %in% Pub_Seurat$hash) # include only cells in the final published meta data
Merged_Liver@meta.data <- left_join(Merged_Liver@meta.data, Pub_Seurat@meta.data, "hash")
rownames(Merged_Liver@meta.data) <- colnames(Merged_Liver)

message("Re-naming meta data columns")
# harmonize column names
Merged_Liver@meta.data <- Merged_Liver@meta.data[ , -which(colnames(Merged_Liver@meta.data) %in% c("orig.ident.x", 
                                                                                                    "patient", 
                                                                                                    "cause_of_disease", 
                                                                                                    "path",
                                                                                                    "Sample",
                                                                                                    "species",
                                                                                                    "nCount_RNA.y",
                                                                                                    "nFeature_RNA.y",
                                                                                                    "nGene",
                                                                                                    "nUMI",
                                                                                                    "percent.mito",
                                                                                                    "condition.y",
                                                                                                    "annotation_original"))]

Merged_Liver$cellType <- Merged_Liver@meta.data$annotation_lineage
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "dataset")] <- "sample"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "aetiology")] <- "cause_of_disease" 
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "liver")] <- "patient"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "annotation_lineage")] <- "cell_Lineage"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "annotation_indepth")] <- "cellType_Specific"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "nCount_RNA.x")] <- "nCount_RNA"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "nFeature_RNA.x")] <- "nFeature_RNA"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "condition.x")] <- "condition"
colnames(Merged_Liver@meta.data)[which(colnames(Merged_Liver@meta.data) == "orig.ident.y")] <- "orig.ident"
Merged_Liver$condition <-  ifelse(Merged_Liver$condition == "healthy", "Uninjured", "Diseased")
Merged_Liver$Pre_Ex <- Merged_Liver$cause_of_disease
Merged_Liver$organ <- "liver"
Merged_Liver$tissue <- "NPC_liver"
Merged_Liver$age <- "DONT_HAVE"
Merged_Liver[["percent.mt"]] <- PercentageFeatureSet(Merged_Liver, pattern = "^MT-")
Merged_Liver[["percent.ribo"]] <- PercentageFeatureSet(Merged_Liver, pattern = "^RP")
Merged_Liver$condition == "fibrotic", "Diseased", "Uninjured"

# more harmonizing




message("Saving merged liver object")
save(Merged_Liver, file = "data/liver/Seurat/Merged_Liver.Rds")

# (4) Now we read in the kidney fibrosis data from Henderson's lab
# Reading in Henderson's kidney scRNA-Seq data
# All samples were FACS sorted
message("Reading in henderson kidney data")
Samples <- list.dirs("data/kidney/matrices", recursive = FALSE, full.names = FALSE)
objectList <- list()
for (projectName in Samples) {
    print(projectName)
    mtx <- readMM(paste0("data/kidney/matrices/", projectName,"/",list.files(paste0("data/kidney/matrices/", projectName))[grep(".mtx", list.files(paste0("data/kidney/matrices/", projectName)))]))

    dim(mtx)
    rowData <- read.table(paste0("data/kidney/matrices/", projectName,"/",list.files(paste0("data/kidney/matrices/", projectName))[grep("rowData", list.files(paste0("data/kidney/matrices/", projectName)))]),
                          sep = ",", header = TRUE)

    rownames(mtx) <- rowData$Gene.Symbol
    dim(rowData)
    colData <- read.table(paste0("data/kidney/matrices/", projectName,"/",list.files(paste0("data/kidney/matrices/", projectName))[grep("colData", list.files(paste0("data/kidney/matrices/", projectName)))]),
                          sep = ",", header = TRUE)
   
    colData <- colData %>% mutate(RN = row_number())
    colData$barcode <- paste(projectName, colData$Patient.ID, colData$RN, sep = "_")
    colnames(mtx) <- colData$barcode
    rownames(colData) <- colData$barcode
    dim(colData)
    objectList[[projectName]] <- CreateSeuratObject(
        counts = mtx, meta.data = colData, min.cells = 0, min.features = 0,
        project = projectName
    )
}

Merged_Kidney <- merge(x = objectList[[1]], y = objectList[c(2:length(objectList))])

# harmonize meta data column names
Merged_Kidney@meta.data <- Merged_Kidney@meta.data[, -which(colnames(Merged_Kidney@meta.data) %in% c("Core.Matrisome.Expression.Score", "RN"))]
colnames(Merged_Kidney@meta.data)[which(colnames(Merged_Kidney@meta.data) == "Annotation.Level.1")] <- "cell_Lineage"
colnames(Merged_Kidney@meta.data)[which(colnames(Merged_Kidney@meta.data) == "Annotation.Level.2")] <- "cellType"
colnames(Merged_Kidney@meta.data)[which(colnames(Merged_Kidney@meta.data) == "Annotation.Level.3")] <- "cellType_Specific"
colnames(Merged_Kidney@meta.data)[which(colnames(Merged_Kidney@meta.data) == "Patient.ID")] <- "patient"
Merged_Kidney$patient <- unlist(lapply(strsplit(Merged_Kidney$patient, split = "_"), getElement, 1))
Merged_Kidney$hash <- Merged_Kidney$barcode

Merged_Kidney$sample <- paste(Merged_Kidney$patient, Merged_Kidney$orig.ident, sep = "_")
Merged_Kidney$organ <- "kidney"
Merged_Kidney$collection <- ifelse(Merged_Kidney$orig.ident == "CD10negative",
    "CD10-",
    ifelse(Merged_Kidney$orig.ident == "Human_CD10plus",
        "CD10+",
        "PDGFRB+"
    )
)
Merged_Kidney$study <- "Kuppe"
Merged_Kidney$tissue <- ifelse(Merged_Kidney$orig.ident == "CD10negative", "non-PT_kidney",
    ifelse(Merged_Kidney$orig.ident == "Human_CD10plus", "PT-kidney",
        "PDGFRB+_sorted_kidney"
    )
)

# include patient level meta data 
Merged_Kidney$condition <- ifelse(Merged_Kidney$patient %in% c(
    "CDm11",
    "CDm5",
    "CDm8",
    "CDm9",
    "CDp1",
    "CDp2",
    "CDp7",
    "Pb2",
    "Pb5",
    "Pb8",
    "Pb9"
),
"Diseased", "Uninjured"
)
Merged_Kidney$cause_of_disease <- ifelse(Merged_Kidney$condition == "Diseased", "CKD", "None")
Merged_Kidney$Pre_Ex <- ifelse(Merged_Kidney$condition == "Diseased", "CKD", "None")
Merged_Kidney$age <- "TBD"
Merged_Kidney$sex <- "TBD"
Merged_Kidney[["percent.mt"]] <- PercentageFeatureSet(Merged_Kidney, pattern = "^MT-")
Merged_Kidney[["percent.ribo"]] <- PercentageFeatureSet(Merged_Kidney, pattern = "^RP")
message("Saving merged kidney Henderson object")
save(Merged_Kidney, file = "data/kidney/Seurat/Merged_Kidney.Rds")

