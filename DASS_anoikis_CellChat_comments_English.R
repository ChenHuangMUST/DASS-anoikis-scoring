library(Seurat)
library(dplyr)
library(patchwork) # Powerful package for combining plots
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
library(future)

load("E:/Anoikis_all/Anoikis_batch_sc_new/Anoikis_batch_cellchat/scRNA_harmony_new.RData")

#### Create a CellChat object ####
cellchat <- createCellChat(object = scRNA_harmony@assays$RNA@data,meta = scRNA_harmony@meta.data,group.by = "celltype")

summary(cellchat)

levels(cellchat@idents)

groupsize <- as.numeric(table(cellchat@idents))

groupsize

#### Load the ligand-receptor database ####
CellChatDB <- CellChatDB.human # Use CellChatDB.mouse if running on mouse data

showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # Use the default CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even when using the whole database

future::plan("multiprocess", workers = 1) # Enable parallel processing

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

#### Compute communication probabilities and infer the cell-cell communication network ####
cellchat <- computeCommunProb(cellchat,population.size = T)

cellchat <- filterCommunication(cellchat, min.cells = 10)## Filter out cell-cell communication networks with fewer than 10 cells

df.net <- subsetCommunication(cellchat) 

#### Infer cell-cell communication at the signaling pathway level ####
cellchat <- computeCommunProbPathway(cellchat)

#### Compute the aggregated cell-cell communication network ####
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight

par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#### Heatmaps of interaction number and strength differences ####
par(mfow=c(1,1))

h1 <- netVisual_heatmap(cellchat)

h2 <- netVisual_heatmap(cellchat,measure = "weight")

h1+h2


#### Pairwise analysis between different groups ####
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)

table(scRNA_harmony$celltype)

table(scRNA_harmony$orig.ident)

Idents(scRNA_harmony) <- "orig.ident"

sco.res <- subset(scRNA_harmony, idents = c("GSM4675273", "GSM4675274", "GSM4675275"))

table(sco.res$celltype)

sco.sen <- subset(scRNA_harmony, idents = c("GSM4675276", "GSM4675277"))

table(sco.sen$celltype)

# Create CellChat objects
cco.res <- createCellChat(sco.res@assays$RNA@data, meta = sco.res@meta.data, group.by = "celltype")

cco.sen <- createCellChat(sco.sen@assays$RNA@data, meta = sco.sen@meta.data, group.by = "celltype")

save(cco.res, cco.sen, file = "cco.rda")

## Cell-cell communication network analysis
# Analyze the cell-cell communication network of the cco.res group
cellchat.res <- cco.res

cellchat.res@DB <- CellChatDB.human

cellchat.res <- subsetData(cellchat.res)

cellchat.res <- identifyOverExpressedGenes(cellchat.res)

cellchat.res <- identifyOverExpressedInteractions(cellchat.res)

cellchat.res <- computeCommunProb(cellchat.res, raw.use = T, population.size = T)

cellchat.res <- computeCommunProbPathway(cellchat.res)

cellchat.res <- aggregateNet(cellchat.res)

cellchat.res <- netAnalysis_computeCentrality(cellchat.res, slot.name = "netP")

cco.res <- cellchat.res

saveRDS(cco.res, file = "cco.res.rds")

# Analyze the cell-cell communication network of the cco.sen group
cellchat.sen <- cco.sen

cellchat.sen@DB <- CellChatDB.human

cellchat.sen <- subsetData(cellchat.sen)

cellchat.sen <- identifyOverExpressedGenes(cellchat.sen)

cellchat.sen <- identifyOverExpressedInteractions(cellchat.sen)

cellchat.sen <- computeCommunProb(cellchat.sen, raw.use = T, population.size = T)

cellchat.sen <- filterCommunication(cellchat.sen)

cellchat.sen <- computeCommunProbPathway(cellchat.sen)

cellchat.sen <- aggregateNet(cellchat.sen)

cellchat.sen <- netAnalysis_computeCentrality(cellchat.sen, slot.name = "netP")

cco.sen <- cellchat.sen

saveRDS(cco.sen, file = "cco.sen.rds")

# Merge CellChat objects
cco.list <- list(res = cco.res, sen = cco.sen)

cellchat.list <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = T)

#### Visualization ####
# Global overview across all cell groups: compare interaction number and strength
gg1 <- compareInteractions(cellchat.list, show.legend = F, group = c(1,2), measure = "count")

gg2 <- compareInteractions(cellchat.list, show.legend = F, group = c(1,2), measure = "weight")

p <- gg1 + gg2

p

# Network plots comparing the number of cell-cell interactions
par(mfrow = c(1,2))

weight.max <- getMaxWeight(cco.list, attribute = c("idents", "count"))

for(i in 1:length(cco.list)){
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge = T, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}

# Network plots comparing the number of interactions among selected cell groups
par(mfrow = c(1,2))

s.cell <- c("T_cell", "Epithelial", "Endothelial","B_cell","Fibroblasts")

count1 <- cco.list[[1]]@net$count[s.cell, s.cell]

count2 <- cco.list[[2]]@net$count[s.cell, s.cell]

weight.max2 <- max(max(count1), max(count2))

netVisual_circle(count1, weight.scale = T, label.edge = T, edge.weight.max = weight.max2, edge.width.max = 12,
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge = T, edge.weight.max = weight.max2, edge.width.max = 12,
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))

## Identify and visualize conserved and group-specific signaling pathways
# Comparative analysis of signaling pathway strength

w1 <- rankNet(cellchat.list, mode = "comparison", stacked = T, do.stat = T)

w2 <- rankNet(cellchat.list, mode = "comparison", stacked = F, do.stat = T)

w <- w1 + w2

w

# Compare cellular signaling patterns
library(ComplexHeatmap)
# Overall signaling pattern comparison
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)

ht1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union,
                                         title = names(cco.list)[1], width = 8, height = 10)

ht2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                         title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Outgoing signaling pattern comparison
ht1_out <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(cco.list)[1], width = 8, height = 10,color.heatmap = "GnBu")

ht2_out <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(cco.list)[2], width = 8, height = 10,color.heatmap = "GnBu")
draw(ht1_out + ht2_out, ht_gap = unit(0.5, "cm"))

# Incoming signaling pattern comparison
ht1_in <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union,
                                         title = names(cco.list)[1], width = 8, height = 10, color.heatmap = "OrRd")

ht2_in <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                         title = names(cco.list)[2], width = 8, height = 10, color.heatmap = "OrRd")

draw(ht1_in + ht2_in, ht_gap = unit(0.5, "cm"))

#### Comparison of specific signaling pathways ####
cco.list$sen@netP$pathways

cco.list$res@netP$pathways

com <- intersect(cco.list$sen@netP$pathways,cco.list$res@netP$pathways)

pathways.show <- c("FN1")

weight.max3 <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)

par(mfrow = c(1,2), xpd = TRUE)

for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max3[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# Heatmap
par(mfrow = c(1,2), xpd = TRUE)

ht_p <- list()

for(i in 1:length(cco.list)){
  ht_p[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ", names(cco.list)[i]))
}

ComplexHeatmap::draw(ht_p[[1]] + ht_p[[2]], ht_gap = unit(0.5, "cm"))

## Ligand-receptor comparison analysis
# Bubble plot showing differences across all ligand-receptor pairs
levels(cellchat.list@idents$joint)

r <- netVisual_bubble(cellchat.list, sources.use = c(2,3), targets.use = c(1,4,8), comparison = c(1,2), angle.x = 45)

r

#### Intersect with differentially expressed genes from the previous six-/three-group analysis ####
load("E:/Anoikis_all/Anoikis_batch_sc_new/Anoikis_batch_cellchat/limma_dif_by_wu2.Rdata")

anno_act_marker <- as.data.frame(deg_list[1])

anno_int_marker <- as.data.frame(deg_list[2])

anno_rep_marker <- as.data.frame(deg_list[3])

anno_act_up <- subset(anno_act_marker,logFC>0 & adj.P.Val< 0.05)

anno_rep_up <- subset(anno_rep_marker,logFC>0 & adj.P.Val< 0.05)

anno_int_up <- subset(anno_int_marker,logFC>0 & adj.P.Val< 0.05)

act_com <- intersect(rownames(anno_act_up),pathway.union)

rep_com <- intersect(rownames(anno_rep_up),pathway.union)

int_com <- intersect(rownames(anno_int_up),pathway.union)

#### APP signaling pathway in the repressed group ####
pathways.show <- c("APP")

weight.max3 <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)

par(mfrow = c(1,2), xpd = TRUE)

for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max3[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

#netAnalysis_contribution(cellchat.res,signaling = pathways.show)

#### Signaling pathway in the intermediate group ####
pathways.show <- c("ADGRE5")

weight.max3 <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)

par(mfrow = c(1,2), xpd = TRUE)

for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max3[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

#### Signaling pathway in the activated group ####
pathways.show <- c("SPP1")

weight.max3 <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)

par(mfrow = c(1,2), xpd = TRUE)

for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max3[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

par(mfrow = c(1,2), xpd = TRUE)

ht_p <- list()

for(i in 1:length(cco.list)){
  ht_p[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                                 title.name = paste(pathways.show, "signaling ", names(cco.list)[i]))
}

ComplexHeatmap::draw(ht_p[[1]] + ht_p[[2]], ht_gap = unit(0.5, "cm"))
