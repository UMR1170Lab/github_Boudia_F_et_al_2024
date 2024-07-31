#' @title annotation_SingleR_dataframe_output
#'
#' @description Cette fonction renvoie un dataframe de meta.data de Seurat contenant une colonne all_single_R_annotation, celle-ci indiquant pour chaque cellule le type cellulaire inféré par SingleR en se basant sur un jeu de données de référence spécifié par l utilisateur
#'
#' @param datatoannot : objet Seurat (pas une intégration) qui peut être employé par SingleR .
#' @param reference : datas qui permettent d'annoter le Single-Cell. Par exemple NovershternHematopoieticData() qui vient du package celldex. En fait, il faut obligatoirement que les données proviennent du package celldex.
#' @return : cette fonction renvoie un dataframe de meta.data Seurat contenant une colonne all_single_R_annotation, celle-ci indiquant pour chaque cellule le type cellulaire inféré par SingleR en se basant sur un jeu de données de référence spécifié par l utilisateur
#' @importFrom SingleR SingleR


annotation_SingleR_dataframe_output <- function(datatoannot,reference = celldex::NovershternHematopoieticData()){
    ## ci-dessus accolade ouvrante de la fonction annotation_SingleR_dataframe_output
    
    require(SingleR)
    
    #annotations de référence
    ref <- reference
    
    
    print(datatoannot)
    
    datatoannotdup <- datatoannot
    print(class(ref))
    #création de l'objet d'annotation SingleR. Pour chaque cellule, une inférence du type cellulaire est réalisée.
    objet_singleR <- SingleR::SingleR(test=datatoannotdup@assays$RNA@counts,
    ref= ref@assays@data$logcounts
    ,labels= ref$label.main,
    de.method="wilcox")
    
    
    # SingleR(test= normalized_data_list[[1]]@assays$RNA@counts,
    #                          ref= refdata@assays@data$logcounts
    #                          ,labels= ref$label.main,
    #                           de.method="wilcox")
    
    #création du vecteur d'annotation qui correspond à une colonne de  datatoannotdup@meta.data$all_single_R_annotation
    datatoannotdup@meta.data$all_single_R_annotation <- objet_singleR@listData$pruned.labels
    #renvoie du dataframe de meta.data contenant la colonne d'annotation SingleR se basant sur le dataset de référence spécifié
    return(datatoannotdup@meta.data)
    
    ## ci-dessous accolade fermante de la fonction annotation_SingleR_dataframe_output
}
