# Chargement des fonctions qui vont être employées pour l annotation des types cellulaires

#' @title association_number_condition
#'
#' @description Cette fonction permet de renvoyer le numéro associé à chaque condition dans les noms de lignes du dataframe de metadata d une intégration effectuée avec Seurat.
#'
#' @param integration Un objet de class Seurat correspondant à une intégration réalisée avec cette méthode. Le dataframe de meta.data doit contenir une colonne orig.ident.
#'
#' @importFrom magrittr %>%

association_number_condition <- function(integration){#accolade ouvrante de la fonction association_number_condition
    
    require(magrittr)
    data.frame(
    
    condition = unique(integration@meta.data$orig.ident),
    number =  (lapply(unique(integration@meta.data$orig.ident), function(y){
        
        rownames(integration@meta.data)[which(integration@meta.data$orig.ident==y)] %>% strsplit("_") %>% lapply(function(z){return(z[[2]])}) %>% unlist() %>% unique()
        
        
        
    }) %>% unlist())
    
    )
    
    
    
}#accolade fermante de la fonction association_number_condition
