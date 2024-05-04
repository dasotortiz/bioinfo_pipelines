#### load libraries & utility function 
library(Seurat)

# inputs
object_path <- snakemake@input[[1]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/scrnaseq_processing_seurat/KOcall_NonTargeting/NORMALIZED_object.rds"

# outputs
umap_res_results_path <- snakemake@output[["umap_res_results"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_results.csv" 

# parameters
assay <- snakemake@params[["assay"]] #"SCT" #"RNA"
metadata <- snakemake@params[["metadata"]] #"condition"
control <- snakemake@params[["control"]] #"untreated"


result_dir <- dirname(umap_res_results_path)
# make directories if not exist
if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

### load data
data <- readRDS(file = file.path(object_path))
DefaultAssay(object = data) <- assay
Idents(object = data) <- metadata


### perform clustering
data2 <- data
all.genes <- rownames(data2)
data2 <- ScaleData(data2, features = all.genes)
data2 <- RunPCA(data2, features = VariableFeatures(object = data2))
data2 <- FindNeighbors(data2, dims = 1:10)
data2 <- FindClusters(data2, resolution = 0.6)


### save objetc with new metadata
saveRDS(data2, file=file.path(umap_res_results_path))
