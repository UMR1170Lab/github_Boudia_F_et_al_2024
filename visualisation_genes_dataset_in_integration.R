visualisation_genes_dataset_in_integration <- function(coordinate_umap_by_datasets,normalized_individual_datasets,gene){#accolade ouvrante de la fonction
  #chargement des pakages nécessaires pour le fonctionnement de la fonction
  require(ggplot2)
  require(purrr)
  
  
  lapply(1:length(normalized_individual_datasets),purrr::possibly(function(z){
    #accolade ouvrante de la fonction du lapply
    
    # on sélectionne les données d'expression dans le jeux de données d'origine z de la list normalized_individual_datasets 
    value_gene_original <- normalized_individual_datasets[[z]]@assays$SCT$data[which(rownames(normalized_individual_datasets[[z]]@assays$SCT$data)== gene),,drop = TRUE]
    
    # on donne des noms aux éléments de ce vecteur.
    names(value_gene_original) <- colnames(normalized_individual_datasets[[z]]@assays$SCT$data)
    
    
    #print(value_gene_original)
    print(any(names(value_gene_original) != coordinate_umap_by_datasets[[z]]$identifiants))
    #on range si nécessaire les valeurs d'expression du gène d'intérêt associé à chaque cellule de telle sorte qu'à chaque cellule soit associé à la fois les bonnes coordonnées de UMAP et la bonne valeur d'expressio pour le gène d'intérêt.
    if(any(names(value_gene_original) != coordinate_umap_by_datasets[[z]]$identifiants)==TRUE){
      
      #instrutions dans le cas où il faut faire un réarrangement
      value_gene_original_correct <- sapply(coordinate_umap_by_datasets[[z]]$identifiants,function(x){return(value_gene_original[which(names(value_gene_original) == x)])})
      
    }else{
      
      value_gene_original_correct <- value_gene_original
      
    }
    
    
    #création du dataframe à partir duquel la figure ggplot2 pourra être généré.
    coordinate_umap_specific <- coordinate_umap_by_datasets[[z]]
    
    #création d'une colonne avec les valeurs d'expression pour le gène d'intérêt
    coordinate_umap_specific[[gene]] <- value_gene_original_correct
    
    #affichage des premières lignes du dataframe qui va permettre la génération du nuage de points ggplot2.
    print(head(coordinate_umap_specific))
    
    #génération de la figure  
    print(
    ggplot(coordinate_umap_specific,aes_string(x = "umap_1",y = "umap_2",color = gene))+
      geom_point(size =0.5)+
      scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral")))+
      ggtitle(paste0(names(normalized_individual_datasets)[z],"_gene_",gene))+
      theme_classic()
    
    #scale_color_viridis_c()
    )   
  },"gene_not_found_or_problems")
  
  
  )
  
  
  
  
  
#accolade fermante de la fonction  
}
