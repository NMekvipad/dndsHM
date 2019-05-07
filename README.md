
# dndsHM

R-script containing function to generate heatmap of pairwise dn/ds of value from dn/ds value calculated by MEGA program and other utilty functions for generating dn/ds heatmap with phylogenetic tree e.g function to generate list of correctly display scientific name (binomial nomenclature with italicized part for species name and non-italicized part for strain name) in plot. 

This script was written for the article titled [Evolution of organellar genes of chlorophyte algae: Relevance to phylogenetic inference](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0216608).


### Prerequisites
MEGA is required for dn and ds value calculation. [MEGA](https://www.megasoftware.net/)

This script depends on:
```
  ape version 5.2, ggdendro version 0.1-20, ggplot2 version 3.1.0, 
  grid version 3.4.4, phylogram version 2.1.0,  reshape2 version 1.4.3,  
  scales version 1.0, stringr version 1.3.1
  
```
### Code example

Lower triangular matrix of dn or ds value from MEGA program can be read and converted into symmetric distance matrix using read_distancematrix function. 

```
dn_sample <- read_distancematrix(dn_file_path)

```
Dendrogram in newick format can be read and converted into dendro type object in R with sorted leaf label using read_dendrogram function. 

```
tree_sample <- read_dendrogram(tree_file_path)

```
The samples of input file can be found in sample_file folder.

MEGA program calculates dn and ds separately, so to generate dn/ds matrix we use dnds_matrix function. The order of label in matrix can be specified by specifying vector of ordered label to the dendrogram_label argument.   

```
dnds_dataframe <- dnds_matrix(dn_mat=dn_sample, ds_mat=ds_sample, dendrogram_label=ordered_label)

```
To generate heatmap as shown in the article, we can use the dndsheatmap function. The function takes 4 arguments. 
    dnds_mat: a symmetric matrix of dn/ds value
    taxon_list:   the list of vector of species in each taxonomic group
    col_palette:  the vector of colors that will be used to color the species name of each     
                  taxonomic group; the order of color in vector are the same as the order of     
                  taxonomic group in taxon_list argument
    ylabels:  the vector of heatmap's y-axis labels

```
# taxon_list example 
taxon_list <- list('taxonomic_group1' = c('species1', 'species2', 'species3', 'species4'),
                   'taxonomic_group2' = c('species5', 'species6', 'species7'),
                   'taxonomic_group3' = c('species8', 'species9', 'species10')
                   )

# col_palette example
## red color corresponding to taxonomic group1
## green color corresponding to taxonomic group2
## blue color corresponding to taxonomic group3

col_palette <- c('red', 'green', 'blue') 

```
To plot dendrogram together with heatmap we use grid.newpage() function 
```
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))
print(dendro_plot, vp = viewport(x = 0.585, y = 0.875, width = 0.524, height = 0.2))
print(heatmap_plot, vp = viewport(x = 0.5, y = 0.42, width = 0.8, height = 1.0))
```

We also provide the italic_sciname function to generate the list of scientific name that can be correctly display in proper text fomatting (binomial nomenclature with italicized part for species name and non-italicized part for strain name). 

The function take the vector of name and the regular expression pattern to regcognize the italicized part and non-italicized part in scientific name.

```
# ([A-Za-z]+?\\s[A-Za-z.]+) to recognize scientific name
# (\\s[SAG|UTEX|UNA|CCMP|UTEX][A-Za-z0-9.\\s]*) to recognize strain name 
pattern <- '([A-Za-z]+?\\s[A-Za-z.]+)(\\s[SAG|UTEX|UNA|CCMP|UTEX][A-Za-z0-9.\\s]*)*'
name_vec <- c("Chlorotetraedron incus SAG43.81", "Neochloris aquatica UTEX B138")
ita_name <- italic_sciname(name_vec, pattern) 

```



### Authors

**Nuttapong Mekvipad** [https://github.com/NMekvipad]

