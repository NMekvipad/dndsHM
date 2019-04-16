library(ape)
library(reshape2)
library(phylogram)
library(ggdendro)
library(ggplot2)
library(scales)
library(grid)
library(stringr)

setwd('/home/nuttapong/Desktop/file_for_run/dnds')

col_palette <- c("deeppink", "darkgreen", "green3", "orangered", "yellow3")
taxon_list <- list('chloro' = c("Bracteacoccus_aerius", "Pseudomuriella_schumacherensis_SAG2137",
                                "Chlorotetraedron_incus_SAG43.81", "Neochloris_aquatica_UTEX_B138",
                                "Chromochloris_zofingiensis_UTEX56"),
                   'ulvo' = c("Gloeotilopsis_planctonica_SAG29.93", "Gloeotilopsis_sarcinoidea_UTEX1710",
                              "Pseudendoclonium_akinetum", "Ulva_fasciata",
                              "Ulva_flexuosa", "Ulva_sp._UNA00071828", "Oltmannsiellopsis_viridis"),
                   'strepto' = c("Mesostigma_viride", "Chaetosphaeridium_globosum",
                                 "Chara_vulgaris", "Chlorokybus_atmophyticus", 
                                 "Entransia_fimbriata_UTEX_LB2353"),
                   'treb' = c("Auxenochlorella_protothecoides", "Chlorella_sp._ArM0029B",
                              "Chlorella_sorokiniana_UTEX1230", "Chlorella_variabilis_NC64A",
                              "Coccomyxa_subellipsoidea_C169"),
                   'pras' = c("Nephroselmis_olivacea", "Ostreococcus_tauri_RCC1561",
                              "Pyramimonas_parkeae_NIES254", "Prasinoderma_coloniale_CCMP1220")
                   )

file_list <- list.files(pattern = 'dn.csv')
gene_list <- str_replace_all(file_list, 'dn.csv', '')
  
for(gene_name in gene_list){
  dn_file <- paste0(gene_name, 'dn.csv', collapse = '')
  dn_dat <- read_distancematrix(dn_file)
  ds_file <- paste0(gene_name, 'ds.csv', collapse = '')
  ds_dat <- read_distancematrix(ds_file)
  tree_file <- paste0(gene_name, '_RAxML_bipartitions.result.newick', collapse = '')
  dendro <- read_dendrogram(tree_file)
  dendro_in <- dendro_data(dendro)
  tree_label <- as.vector(dendro_in$labels$label)
  
  dnds_df <- dnds_matrix(dn_dat, ds_dat, tree_label)
  dnds_mat <- as.matrix(dnds_df)
  
  ylabels_preformatted <- rev(str_replace_all(colnames(dnds_mat), '_', ' '))
  pattern <- '([A-Za-z]+?\\s[A-Za-z.]+)(\\s[SAG|UTEX|UNA|CCMP|UTEX|ArM|NC|C|RCC|NIES][A-Za-z0-9.\\s]*)*'
  ylabels <- italic_sciname(ylabels_preformatted, pattern)
  heatmap_plot <- dndsheatmap(dnds_mat, col_palette, taxon_list, ylabels)
  dendro_plot <- plotdendro(dendro_in)
  
  file_name = paste0(gene_name, '.pdf', collapse = '')
  # all value here need to be fixed to get perfectly align plot
  pdf(file_name, width = 11, height = 7)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2,1)))
  print(dendro_plot, vp = viewport(x = 0.585, y = 0.875, width = 0.524, height = 0.2))
  print(heatmap_plot, vp = viewport(x = 0.5, y = 0.42, width = 0.8, height = 1.0))
  dev.off()
}

# cox1 need to replot to adjust width of dendrogram to width = 0.533 and x to 0.59

gene_name <- 'cox1'

dn_file <- paste0(gene_name, 'dn.csv', collapse = '')
dn_dat <- read_distancematrix(dn_file)
ds_file <- paste0(gene_name, 'ds.csv', collapse = '')
ds_dat <- read_distancematrix(ds_file)
tree_file <- paste0(gene_name, '_RAxML_bipartitions.result.newick', collapse = '')
dendro <- read_dendrogram(tree_file)
dendro_in <- dendro_data(dendro)
tree_label <- as.vector(dendro_in$labels$label)

dnds_df <- dnds_matrix(dn_dat, ds_dat, tree_label)
dnds_mat <- as.matrix(dnds_df)

ylabels_preformatted <- rev(str_replace_all(colnames(dnds_mat), '_', ' '))
pattern <- '([A-Za-z]+?\\s[A-Za-z.]+)(\\s[SAG|UTEX|UNA|CCMP|UTEX|ArM|NC|C|RCC|NIES][A-Za-z0-9.\\s]*)*'
ylabels <- italic_sciname(ylabels_preformatted, pattern)
heatmap_plot <- dndsheatmap(dnds_mat, col_palette, taxon_list, ylabels)
dendro_plot <- plotdendro(dendro_in)

file_name = paste0(gene_name, '.pdf', collapse = '')
# all value here need to be fixed to get perfectly align plot
pdf(file_name, width = 11, height = 7)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))
print(dendro_plot, vp = viewport(x = 0.59, y = 0.875, width = 0.533, height = 0.2))
print(heatmap_plot, vp = viewport(x = 0.5, y = 0.42, width = 0.8, height = 1.0))
dev.off()


# the rest in following list need to replot to adjust width of dendrogram to width = 0.517 and 
# x=0.583

gene_list <- c('cox2', 'nad3', 'nad4L', 'nad5', 'petB', 'psaA', 'psaB', 'psbC', 'psbD',
               'rpl5', 'rpl20', 'rps14')


for(gene_name in gene_list){
  dn_file <- paste0(gene_name, 'dn.csv', collapse = '')
  dn_dat <- read_distancematrix(dn_file)
  ds_file <- paste0(gene_name, 'ds.csv', collapse = '')
  ds_dat <- read_distancematrix(ds_file)
  tree_file <- paste0(gene_name, '_RAxML_bipartitions.result.newick', collapse = '')
  dendro <- read_dendrogram(tree_file)
  dendro_in <- dendro_data(dendro)
  tree_label <- as.vector(dendro_in$labels$label)
  
  dnds_df <- dnds_matrix(dn_dat, ds_dat, tree_label)
  dnds_mat <- as.matrix(dnds_df)
  
  ylabels_preformatted <- rev(str_replace_all(colnames(dnds_mat), '_', ' '))
  pattern <- '([A-Za-z]+?\\s[A-Za-z.]+)(\\s[SAG|UTEX|UNA|CCMP|UTEX|ArM|NC|C|RCC|NIES][A-Za-z0-9.\\s]*)*'
  ylabels <- italic_sciname(ylabels_preformatted, pattern)
  heatmap_plot <- dndsheatmap(dnds_mat, col_palette, taxon_list, ylabels)
  dendro_plot <- plotdendro(dendro_in)
  
  file_name = paste0(gene_name, '.pdf', collapse = '')
  # all value here need to be fixed to get perfectly align plot
  pdf(file_name, width = 11, height = 7)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2,1)))
  print(dendro_plot, vp = viewport(x = 0.583, y = 0.875, width = 0.517, height = 0.2))
  print(heatmap_plot, vp = viewport(x = 0.5, y = 0.42, width = 0.8, height = 1.0))
  dev.off()
}


