library(liana)
library(nichenetr)
library(Seurat)
library(ggrepel)
library(cowplot)
library(dplyr)

ligand_target_matrix <- readRDS("ligand_target_matrix.rds")

clist = c('CID3941')#,'CID44991','CID4513','CID4515','CID4523','CID4530N','CID4535')#'CID3941','CID3946','CID3948','CID3963','CID4040','CID4067','CID4290A','CID4398','CID44041','CID4461','CID4463','CID4465','CID4471','CID4495','CID44971','CID44991','CID4513','CID4515','CID4523','CID4530N','CID4535')
for (cn in clist) {
  data <- Read10X(data.dir = paste("./patient_data/", cn ,"/", sep = ""),gene.column=1)
  metadata <- read.csv(paste("./patient_data/", cn ,"/metadata.csv", sep = ""))
  rownames(metadata) = metadata$X
  expression = data@x
  geneset_oi = readr::read_tsv("mesenchymal.txt", col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)]
  background_expressed_genes = data@Dimnames[1][[1]] %>% .[. %in% rownames(ligand_target_matrix)]
  potential_ligands = readr::read_tsv("input_nichenet_genelist_lum_lig.txt", col_names = "gene") %>% pull(gene) %>% .[. %in% colnames(ligand_target_matrix)]
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
}

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
