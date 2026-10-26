library(ggrepel)
library(pheatmap)
library(tidydr)
library(Seurat)
library(celldex)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(SingleR)
library(clustree)
library(cowplot)
library(SCpubr)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)

library(ggpubr)
library(ggplot2)
library(ggalluvial)
library(ComplexHeatmap)
library(scCustomize)
library(scTenifoldKnk)


file_paths <- c(
  "GSM8245026_ctrl_05_2_barcodes.tsv.gz",
  "GSM8245026_ctrl_05_2_features.tsv.gz",
  "GSM8245026_ctrl_05_2_matrix.mtx.gz",
  "GSM8245017_ctrl_01_1_barcodes.tsv.gz",
  "GSM8245017_ctrl_01_1_features.tsv.gz",
  "GSM8245017_ctrl_01_1_matrix.mtx.gz",
  "GSM8245018_ctrl_01_2_barcodes.tsv.gz",
  "GSM8245018_ctrl_01_2_features.tsv.gz",
  "GSM8245018_ctrl_01_2_matrix.mtx.gz",
  "GSM8245019_ctrl_02_1_barcodes.tsv.gz",
  "GSM8245019_ctrl_02_1_features.tsv.gz",
  "GSM8245019_ctrl_02_1_matrix.mtx.gz",
  "GSM8245020_ctrl_02_2_barcodes.tsv.gz",
  "GSM8245020_ctrl_02_2_features.tsv.gz",
  "GSM8245020_ctrl_02_2_matrix.mtx.gz",
  "GSM8245021_ctrl_03_1_barcodes.tsv.gz",
  "GSM8245021_ctrl_03_1_features.tsv.gz",
  "GSM8245021_ctrl_03_1_matrix.mtx.gz",
  "GSM8245022_ctrl_03_2_barcodes.tsv.gz",
  "GSM8245022_ctrl_03_2_features.tsv.gz",
  "GSM8245022_ctrl_03_2_matrix.mtx.gz",
  "GSM8245023_ctrl_04_1_barcodes.tsv.gz",
  "GSM8245023_ctrl_04_1_features.tsv.gz",
  "GSM8245023_ctrl_04_1_matrix.mtx.gz",
  "GSM8245024_ctrl_04_2_barcodes.tsv.gz",
  "GSM8245024_ctrl_04_2_features.tsv.gz",
  "GSM8245024_ctrl_04_2_matrix.mtx.gz",
  "GSM8245025_ctrl_05_1_barcodes.tsv.gz",
  "GSM8245025_ctrl_05_1_features.tsv.gz",
  "GSM8245025_ctrl_05_1_matrix.mtx.gz"
)


sample_names <- basename(file_paths)

seurat_objects <- list()

for (i in seq_along(file_paths)) {
  file_path <- file_paths[i]
  
  
  mat <- Read10X(data.dir = file_path)
  
  
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200,
    names.field = NULL,
    names.delim = NULL
  )
  
  
  seurat_objects[[i]] <- seurat_obj
}


mmRNA <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = sample_names
)


table(mmRNA@meta.data$orig.ident)


save_path <- "mmRNA.rds"
saveRDS(mmRNA, file = save_path)


mmRNA <- readRDS("mmRNA.rds")








mmRNA[["percent.mt"]] <- PercentageFeatureSet(mmRNA, pattern = "^MT-")


HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB.genes <- intersect(HB.genes, rownames(mmRNA))
mmRNA[["percent.HB"]] <- PercentageFeatureSet(mmRNA, features = HB.genes)

View(mmRNA@meta.data)


p <- VlnPlot(mmRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4)
ggsave("vln_plot_before_human.pdf", plot = p, width = 10, height = 8)




mmRNA <- subset(mmRNA, subset =
                  nCount_RNA >= 500 & nCount_RNA <= 50000 &
                  nFeature_RNA >= 250 & nFeature_RNA <= 5000 &
                  percent.mt <= 15 &
                  percent.HB <= 1
)

VlnPlot(mmRNA , features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
ggsave("vln_plot_after.pdf", plot = last_plot(), width = 10, height = 8)
dev.off()





mmRNA <- NormalizeData(mmRNA, normalization.method = "LogNormalize", scale.factor = 10000)
mmRNA <- FindVariableFeatures(mmRNA, selection.method = "vst", nfeatures = 2000)
mmRNA <- ScaleData(mmRNA, features = VariableFeatures(mmRNA))



mmRNA <- RunPCA(mmRNA, npcs = 50, verbose = FALSE)


c1 <- DimPlot(mmRNA, reduction = "pca", group.by = "orig.ident") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    legend.position = "right"
  )

c1

ggsave("pca_nonharmony.pdf", plot = c1, width = 7.5, height = 6)



ElbowPlot(mmRNA, ndims = 5)


chosen_dims <- 1:5
mmRNA_harmony <- RunHarmony(
  object = mmRNA,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = chosen_dims,
  verbose = FALSE
)
c2 <- DimPlot(mmRNA_harmony, reduction = "harmony", group.by = "orig.ident") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
    legend.position = "right"
  )

c2



mmRNA_harmony <- FindNeighbors(mmRNA_harmony, reduction = "harmony", dims = chosen_dims)


for (res in c(0.1, 0.2, 0.3, 0.5, 0.8)) {
  mmRNA_harmony <- FindClusters(mmRNA_harmony, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}


p_tree <- clustree(mmRNA_harmony@meta.data, prefix = "RNA_snn_res.")
p_tree
ggsave("cluster_tree_大群.pdf", plot = p_tree, width = 10, height = 10)


mmRNA_harmony <- FindClusters(mmRNA_harmony, resolution = 0.5)

mmRNA_harmony <- RunUMAP(mmRNA_harmony, reduction = "harmony", dims = chosen_dims)



pdf(file = "After_cluster_UMAP.pdf", width = 7, height = 6)
DimPlot(mmRNA_harmony, label = TRUE, repel = TRUE,
        label.size = 4, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = "right") +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

dev.off()






mmRNA_harmony.markers <- FindAllMarkers(
  mmRNA_harmony,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  group.by = "seurat_clusters"
)
head(mmRNA_harmony.markers)
mmRNA_harmony.markers$pct.3 <- abs(mmRNA_harmony.markers$pct.1 - mmRNA_harmony.markers$pct.2)
head(mmRNA_harmony.markers)


top10 <- mmRNA_harmony.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = pct.3)
write.csv(top10, file = "top10FindAllMarkers.csv")
write.csv(mmRNA_harmony.markers, file = "FindAllMarkers.csv")



mmRNA_harmony$Celltype <- "Unknown"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(0, 9)] <- "T cells"



mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(1, 3, 6, 7)] <- "Macrophages"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(2, 11, 16)] <- "B cells"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(4, 8, 12, 15)] <- "Cancer Cells"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(5)] <- "Myofibroblasts"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(10)] <- "Fibroblasts"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(13)] <- "ECs"


mmRNA_harmony$Celltype[mmRNA_harmony$seurat_clusters %in% c(14)] <- "Plasma Cells"



table(mmRNA_harmony$Celltype, useNA = "ifany")


table(mmRNA_harmony$seurat_clusters, mmRNA_harmony$Celltype)




cluster_colors <- c("#EDCFD8", "#CBE0A7", "#c0504d", "#C98F96", "#9AD1C3",
                    "#4f80bd", "#EFB27C", "#696969", "#67A59B", "#C0937E",
                    "#D9D9D9", "#F4B400")
p <- DimPlot(mmRNA_harmony,
             reduction = "umap",
             group.by = "Celltype",
             label = TRUE,
             label.size = 4,
             cols = cluster_colors,
             repel = TRUE) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())


p



table(mmRNA_harmony@meta.data$Celltype)

Hu_SMC <- mmRNA_harmony[,mmRNA_harmony$Celltype %in% c('Macrophages')]

DefaultAssay(Hu_SMC) <- "RNA"
Hu_SMC <- NormalizeData(Hu_SMC, normalization.method = "LogNormalize", scale.factor = 1e4)
Hu_SMC <- FindVariableFeatures(Hu_SMC, selection.method = "vst", nfeatures = 2000)
Hu_SMC <- ScaleData(Hu_SMC, features = VariableFeatures(Hu_SMC))

Hu_SMC <- RunPCA(Hu_SMC, features = VariableFeatures(Hu_SMC), npcs = 30, verbose = FALSE)

ElbowPlot(Hu_SMC, ndims = 30)


chosen_dims <- 1:7
Hu_SMC <- RunHarmony(
  object = Hu_SMC,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = chosen_dims,
  verbose = FALSE
)

Hu_SMC <- FindNeighbors(Hu_SMC, reduction = "harmony", dims = chosen_dims)

for (res in c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8)) {
  Hu_SMC <- FindClusters(Hu_SMC, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}


p_tree <- clustree(Hu_SMC@meta.data, prefix = "RNA_snn_res.")
p_tree
ggsave("cluster_tree_巨噬.pdf", plot = p_tree, width = 10, height = 10)


Hu_SMC <- FindClusters(Hu_SMC, resolution = 0.2)



Hu_SMC <- RunUMAP(Hu_SMC, reduction = "harmony", dims = chosen_dims)

cluster_colors <- c("#EDCFD8", "#CBE0A7", "#c0504d", "#9AD1C3", "#C98F96",
                    "#4f80bd", "#EFB27C", "#696969", "#67A59B", "#C0937E",
                    "#D9D9D9", "#F4B400")

p <- DimPlot(Hu_SMC, label = TRUE, repel = TRUE,
             group.by = "seurat_clusters",
             label.size = 5, pt.size = 0.5,
             cols = cluster_colors
) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = "right") +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
p
ggsave("umap_注释巨噬细胞.pdf", plot = p, width = 7, height = 6)




table(Hu_SMC$group)
table(Hu_SMC$orig.ident)
Idents(Hu_SMC) <- "group"

deg_markers <- FindMarkers(
  object = Hu_SMC,
  ident.1 = "BC",
  ident.2 = "Control",
  only.pos = FALSE,
  test.use = "t",
  min.pct = 0,
  logfc.threshold = 0
)
head(deg_markers)
deg_markers$pct.3 <- abs(deg_markers$pct.1 - deg_markers$pct.2)
head(deg_markers)

deg_markers <- deg_markers[deg_markers$p_val_adj < 0.05, ]
deg_markers <- deg_markers[complete.cases(deg_markers), ]

head(deg_markers)


pal <- viridis(n = 10, option = "D", direction = -1)
p <- FeaturePlot_scCustom(
  seurat_object = mmRNA_harmony,
  features = c("SRC"),
  reduction = "umap",
  split.by = "group",
  na_color = "lightgrey",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.8,
  
) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

p







library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggplot2)


DimPlot(Hu_SMC,label = T)


set.seed(12345)
table(Hu_SMC$seurat_clusters)


Hu_SMC$Celltype <- "Unknown"


Hu_SMC$Celltype[Hu_SMC$seurat_clusters %in% c(0)] <- "Macrophages_0"
Hu_SMC$Celltype[Hu_SMC$seurat_clusters %in% c(1)] <- "Macrophages_1"
Hu_SMC$Celltype[Hu_SMC$seurat_clusters %in% c(2)] <- "Macrophages_2"
Hu_SMC$Celltype[Hu_SMC$seurat_clusters %in% c(3)] <- "Macrophages_3"
Hu_SMC$Celltype[Hu_SMC$seurat_clusters %in% c(3)] <- "Macrophages_4"
Hu_SMC$Celltype[Hu_SMC$seurat_clusters %in% c(3)] <- "Macrophages_5"

Idents(Hu_SMC) <- "Celltype"
seu_obj = subset(Hu_SMC,idents=c("Macrophages_0","Macrophages_1","Macrophages_2","Macrophages_3","Macrophages_4","Macrophages_5"))
table(Hu_SMC$Celltype)


data <- GetAssayData(seu_obj, assay = 'RNA', slot = 'counts')
cell_metadata <- seu_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)


cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(cds, preprocess_method = "PCA")


p1=plot_cells(cds, color_cells_by="Celltype",cell_size=0.5,group_label_size=5)
p1

plot_cells(cds, genes=c("CRABP2"),cell_size=0.5,group_label_size=5)



cds<-cluster_cells(cds,resolution=1e-5)
plot_cells(cds)
plot_cells(cds,color_cells_by="partition",
           group_cells_by="partition")



cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Hu_SMC, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed



cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = 'louvain')
cds <- learn_graph(cds)

cds <- learn_graph(cds,
                   verbose=T,
                   learn_graph_control=list(minimal_branch_len=20
                   ))


saveRDS(cds, "巨噬细胞的拟时序monocle3.rds")
cds <- readRDS("巨噬细胞的拟时序monocle3.rds")

save.image(file = "巨噬细胞的拟时序monocle3.RData")
load("拟时序monocle3小胶质.RData")


p = plot_cells(cds, color_cells_by = "Celltype", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
p
ggsave("Monocle3_Epi_Cancer.pdf", p, width = 5, height = 4)

plot_cells(cds,
           color_cells_by = 'Celltype',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,
           group_label_size=4) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())


cluster_colors_zjh <- c("#CBE0A7", "#C0937E", "#c7e0ed", "#EFB27C", "#f05064",
                        "#f05289ff", "#9AD1C3")
cluster_colors <- c("#EDCFD8", "#CBE0A7", "#c0504d", "#C98F96", "#9AD1C3",
                    "#4f80bd", "#EFB27C", "#696969", "#67A59B", "#C0937E")
cluster_colors_FHX <- c("#CBE0A7", "#c7e0ed", "#c0504d", "#EFB27C", "#4f80bd",
                        "#9AD1C3", "#67A59B")
c("Microglia_0","Microglia_1","Microglia_2","Microglia_3","Microglia_4")
custom_colors <- c("Microglia_0" = "#CBE0A7", "Microglia_1" = "#c7e0ed", "Microglia_2" = "#c0504d", "Microglia_3" = "#EFB27C", "Microglia_4" = "#4f80bd")
plot_cells(cds,
           color_cells_by = "Celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=4) +
  scale_color_manual(values = custom_colors) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())


cds <- order_cells(cds)

plot_cells(cds, label_cell_groups = F,
           color_cells_by = "pseudotime",
           label_leaves = F,
           label_branch_points = F,
           graph_label_size = 4,
           cell_size=0.5,
           trajectory_graph_segment_size = 2) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

get_earliest_principal_node <- function(cds, my_select="Microglia_1"){
  cell_ids <- which(colData(cds)[, "Celltype"] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))



plot_cells(cds, label_cell_groups = F,
           color_cells_by = "pseudotime",
           label_leaves = F,
           label_branch_points = F,
           graph_label_size = 4,
           cell_size=0.5,
           trajectory_graph_segment_size = 2) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  scale_color_gradient2(
    low = "grey",
    mid = "#c7e0ed",
    high = "#f05289ff",
    midpoint = 5.5,
    
    name = "Expression"
  ) +
  theme(legend.position = "right")
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

saveRDS(cds, "巨噬细胞的拟时序 monocle3.rds")
cds <- readRDS("Microglia 拟时序 monocle3.rds")




cds_gene <- graph_test(cds,
                       neighbor_graph="principal_graph", cores=4)

write.csv(cds_gene, file = "cds_gene.csv")

cds_gene=read.csv("cds_gene.csv",header = T,row.names = 1)

res_ids <- row.names(subset(cds_gene, q_value < 0.01))
gene_module_df <- find_gene_modules(cds[res_ids,],
                                    resolution=c(10^seq(-6,-1)))
write.csv(gene_module_df, file = "gene_module_df.csv")
gene_module_df=read.csv("gene_module_df.csv",header = T,row.names = 1)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),
                                cell_group=colData(cds)$Celltype)
agg_mat <- aggregate_gene_expression(cds,
                                     gene_module_df,
                                     cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))

pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")



p1 <- p1 +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1.0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  labs(x = "Pseudotime", y = "Expression")

print(p1)


p2 <- plot_genes_in_pseudotime(
  cds[my_genes, ],
  color_cells_by = "pseudotime",
  min_expr = 0.5,
  ncol = 2
)



p2 + scale_color_gradient2(        low = "grey",
                                   mid = "#c7e0ed",
                                   high = "#f05289ff",
                                   name = "Pseudotime") +theme_minimal() +
  theme(
    
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),) +
  
  labs(
    x = "pseudotime",
    y = "Expression",
  )


p2_styled <- p2 +
  
  scale_color_viridis(option = "plasma", name = "Pseudotime") +
  
  theme(
    
    axis.line = element_line(color = "black", linewidth = 2),
    
    axis.ticks = element_line(color = "black", linewidth = 1),
    
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "pseudotime",
    y = "Expression"
  )

p2_styled

my_colors <- c("#440154", "#413985", "#30688D", "#1F928B", "#35B779", "#91D142", "#F8E620")
p2 + scale_color_gradientn(colors = my_colors, name = "Pseudotime")