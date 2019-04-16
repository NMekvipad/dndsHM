library(ape)
library(reshape2)
library(phylogram)
library(ggdendro)
library(ggplot2)
library(scales)
library(stringr)

expand_matrix <- function(data, row_name, col_name = row_name){
  out_df <- data.frame(matrix(nrow = length(row_name), ncol = length(col_name)))
  row.names(out_df) <- row_name
  colnames(out_df) <- col_name
  
  for(i_lab in row_name){
    for(j_lab in col_name){
      val <- data[data[,1] == i_lab & data[,2] == j_lab, 3]
      
      if(length(val) == 0){
        val <- NA
      }
      out_df[i_lab, j_lab] <- val
      
    }
    return(out_df)
  }
  
}

create_symmetric <- function(data){
  data[upper.tri(data)] <- t(data)[upper.tri(data)]
  return(data)
}

read_distancematrix <- function(file){
  dat <- read.table(file, sep = ',', row.names = 1, na.strings=c("","NA","?"))
  symm_dat <- create_symmetric(dat)
  colnames(symm_dat) <- row.names(symm_dat)
  return(symm_dat)
}

read_dendrogram <- function(tree_file){
  
  tree <- read.tree(file = tree_file)
  
  # sort tree by species name to generate good looking tree
  label <- tree$tip.label
  sorted_label <- sort(label)
  rotated_tree <- rotateConstr(tree, sorted_label)
  
  # convert to base r dendrogram
  dendro_sorted <- as.cladogram(as.dendrogram.phylo(rotated_tree))

  return(dendro_sorted)
}

dnds_matrix <- function(dn_mat, ds_mat, dendrogram_label){
  
  dnds_mat <- dn_mat/ds_mat
  # sort row and column
  dnds_mat <- dnds_mat[dendrogram_label, ]
  dnds_mat <- dnds_mat[ ,dendrogram_label]
  dnds_mat[is.na(dnds_mat)] <- 0
  return(dnds_mat)
}

colvec_gen <- function(label, cols, taxon_list){
  cols_vec <- rep('black', length(label))
  for(i in 1:length(taxon_list)){
    bool_vec <- label %in% taxon_list[[i]]
    cols_vec[bool_vec] <- cols[i]
  }
  return(cols_vec)
}


dndsheatmap <- function(dnds_mat, col_palette, taxon_list, ylabels){
  min_dnds <- min(dnds_mat)
  max_dnds <- max(dnds_mat)
  melt_mat <- melt(dnds_mat) # Var1 and Var2 are species name of each pair and value is dn/ds
  col_lab <- colvec_gen(melt_mat$Var1, col_palette, taxon_list)
  col_lab <- rev(col_lab)
  
  if(max_dnds<1){
    heatmap_plot <- ggplot(data = melt_mat, aes(Var1,Var2)) + 
      geom_tile(aes(fill=value)) + scale_y_discrete(name="", limits = rev(levels(melt_mat$Var2)),
                                                    label = ylabels) + 
      scale_fill_gradientn(colours = c("black","red","white"), values = rescale(c(0, max_dnds, 1))) + 
      theme(axis.text.y = element_text(colour = col_lab), 
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
      coord_equal() + 
      guides(fill = guide_colorbar(title = "dn/ds", nbin = 50, 
                                   raster = TRUE, draw.ulim = TRUE, draw.llim = TRUE))
  } else{
    heatmap_plot <- ggplot(data = melt_mat, aes(Var1,Var2)) + 
      geom_tile(aes(fill=value)) + scale_y_discrete(name="", limits = rev(levels(melt_mat$Var2)),
                                                    label = ylabels) + 
      scale_fill_gradientn(colours = c("black","red","white","blue"), values = rescale(c(0, min_dnds, 1, max_dnds))) + 
      theme(axis.text.y = element_text(colour = col_lab), 
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
      coord_equal() + 
      guides(fill = guide_colorbar(title = "dn/ds", nbin = 50, 
                                   raster = TRUE, draw.ulim = TRUE, draw.llim = TRUE))
  }
  return(heatmap_plot)
}

plotdendro <- function(dendrogram_info){
  # convert dendro gram to df for ggplot
  segment_data <- with(
    segment(dendrogram_info), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  
  dendro_plot <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) 
  return(dendro_plot)
}


italic_sciname <- function(name_vec, pattern){
  match_df <- str_match(name_vec, pattern)
  formatted_name <- c()
  for(i in 1:nrow(match_df)){
    if(is.na(match_df[i, 3])){
      formatted_str <- bquote(paste(italic(.(match_df[i, 2]))))
      formatted_name <- c(formatted_name, formatted_str)
    }
    else{
      formatted_str <- bquote(paste(italic(.(match_df[i, 2])), .(match_df[i, 3])))
      formatted_name <- c(formatted_name, formatted_str)
    }
  }
  return(formatted_name)
}
