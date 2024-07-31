#' @title plotclusterpoint_in_original_dataset

#' @description After an integration of several single-cell datasets by Seurat, this function allow to see where the droplets of a cluster of an integration of datasets are in the original datasets

#' @param integrated_data objet Seurat object resulting from the integration of several datasets.

#' @param normalized_data : list of the normalized Seurat object associated to the unnormalised Seurat object that have been used to perform the integration of datasets. Each element of the list should be #' named.
#'
#' @param reduction_normalized_data : list of element of class character indicating the name of the reduction (generally under the form kind of normalisation_technic of dimension reduction_number of dimension kept_representation that is used then). for example SCT_pca_41_umap". The name of list element should be the same than nthe list normalized_data.
#'
#' @export
#' @importFrom Seurat DimPlot
#' @importFrom purrr pmap
#' @importFrom magrittr %>%



plotclusterpoint_in_original_dataset <- function(integrated_data,normalized_data,reduction_normalized_data){
    
    #chargement des packages nécessaires à l utilisation de la fonction
    require(purrr)
    require(Seurat)
    require(magrittr)
    
    ################# Génération d'une liste composé de listes associées à chacuns des clusters déterminés par Seurat.###############
    compclusterdroplets <- lapply(levels(integrated_data@meta.data$seurat_clusters),FUN=function(x){
        
        print(x)
        datatoret <- lapply(as.factor(levels(as.factor(integrated_data@meta.data$orig.ident))),FUN=function(y){
            
            toreturn <- rownames(integrated_data@meta.data)[which(integrated_data@meta.data$orig.ident==y & integrated_data@meta.data$seurat_clusters==x)]  %>% strsplit("_") %>% lapply(FUN=function(x){return(x[1])}) %>% unlist()
            
            return(toreturn)
            
            
        })
        
        names(datatoret) <- as.character(as.factor(levels(as.factor(integrated_data@meta.data$orig.ident))))
        
        
        return(datatoret)
    })
    # On associe un nom à chaque élément de la list générée ci-dessus
    names(compclusterdroplets) <- paste0("cluster_",levels(integrated_data@meta.data$seurat_clusters))
    
    #Récupération des noms des goutelettes qui sont associées à chaque dataset de départ.
    
    dropincluster <- lapply(levels(as.factor(integrated_data@meta.data$orig.ident)),FUN = function(x){
        
        
        
        return(lapply(compclusterdroplets,function(z){
            return(z[[x]])
            
        }))
        
        names(dropincluster) <- levels(as.factor(integrated_data@meta.data$orig.ident))
        
        
    })
    names(dropincluster) <- levels(as.factor(integrated_data@meta.data$orig.ident))
    for(i in 1:nlevels(as.factor(integrated_data@meta.data$orig.ident))){
        
        purrr::pmap(.l=list(data = dropincluster[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]],
        nom=as.list(names(dropincluster[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]]))),
        .f=function(data,nom){
            
            # print(levels(as.factor(integrated_data@meta.data$orig.ident))[i])
            
            normalized_data[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]]@meta.data[[nom]] <- ((rownames(normalized_data[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]]@meta.data) %in% data) %>% as.factor())
            
            print(levels(as.factor(integrated_data@meta.data$orig.ident)))
            print(levels(as.factor(integrated_data@meta.data$orig.ident))[i])
            print(normalized_data[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]]@reductions)
            print(head(normalized_data[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]]@meta.data))
            print(names(normalized_data))
            print(Seurat::DimPlot(normalized_data[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]],
            reduction = reduction_normalized_data[[levels(as.factor(integrated_data@meta.data$orig.ident))[i]]],
            group.by = nom,
            split.by = "orig.ident"))
            
            
            
            
        })
    }
    
}





