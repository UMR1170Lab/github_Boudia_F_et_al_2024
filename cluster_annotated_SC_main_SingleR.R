#' @title cluster_annotated_SC_main_SingleR
#'
#' @description Cette fonction renvoie une représentation du dataset d'entrée en UMAP avec indication des clusters et des annotations  (type cellulaire) qui leurs sont
#' associés (annotation de type cellulaire assez basiques). Elle renvoie également une heatmap représentant les probabilité associées à chaque type cellulaire pour chaque goutelette. Les packages Seurat #'  et SingleR sont employés dans cette fonction.
#'
#' @param datatoannot : objet Seurat qui pourra être utilisé pour les analyses.
#' @param reference : datas qui permettent d'annoter le Single-Cell. Par exemple NovershternHematopoieticData() qui vient du package celldex. En fait, il faut obligatoirement que les données proviennent du package celldex.
#' @param vec_color : vecteur des couleurs qui sont associées à chaque type cllulaire possible de la référence (label.main). par défaut couleurs associées à NovershternHematopoieticData()$label.main

#' @title_rep : chaîne de caractère qui donne le titre que l'on veut donner à la figure pour une intégration donnée
#' @return : cette fonction ne retourne rien en soit. Cette fonction renvoie une représentation du dataset d'entrée en UMAP avec indication des clusters et des annotations  (type cellulaire) qui leurs
#'   sont associés (annotation de type cellulaire assez basiques). Elle renvoie également une heatmap représentant les probabilité associées à chaque type cellulaire pour chaque goutelette. Les packages  #'   Seurat et SingleR sont employés dans cette fonction.

#' @importFrom SingleR SingleR
#' @importFrom Seurat DimPlot
#' @importFrom Seurat DoHeatmap
#' @importFrom ggplot2 theme


cluster_annotated_SC_main_SingleR <- function(datatoannot,reference = celldex::NovershternHematopoieticData(),vec_color = c("Basophils" = "palegreen" , "B cells" = "gold" , "CMPs" = "darkseagreen4", "Dendritic cells" = "purple" , "Eosinophils" = "sienna1" , "Erythroid cells" = "red" , "GMPs" = "grey","Granulocytes" = "orange","HSCs" = "black", "Megakaryocytes" = "blue" , "MEPs" = "pink", "Monocytes" = "aquamarine4", "NK cells" = "antiquewhite3", "NK T cells" = "bisque4",      "CD8+ T cells" = "green",   "CD4+ T cells" = "navajowhite3" ),title_rep){#accolade ouvrante de la fonction cluster_annotated_SC_main
    
    #annotations de référence
    ref <- reference
    
    
    print(datatoannot)
    
    datatoannotdup <- datatoannot
    
    #création de l'objet d'annotation singleR
    objet_singleR <- SingleR::SingleR(test=datatoannotdup@assays$RNA@counts,
    ref= ref@assays@data$logcounts
    ,labels= ref$label.main,
    de.method="wilcox")
    
    
    # SingleR(test= normalized_data_list[[1]]@assays$RNA@counts,
    #                          ref= refdata@assays@data$logcounts
    #                          ,labels= refdata$label.main,
    #                           de.method="wilcox")
    
    #création du vecteur d'annotation
    datatoannotdup@meta.data$all_single_R_annotation <- objet_singleR@listData$pruned.labels
    
    print(
    #Génération de la représentation graphique
    Seurat::DimPlot(datatoannotdup,group.by = "all_single_R_annotation",pt.size =1.4,cols = vec_color)+
    ggplot2::theme(legend.text = element_text(size = 10))
    )
    
    print(
    Seurat::DimPlot(datatoannotdup,
    pt.size = 1.4,label=TRUE,label.size = 8,label.color = "red")+
    ggplot2::theme(legend.text = element_text(size=10))
    
    )
    
    matrice_score <- objet_singleR@listData$scores
    
    rownames(matrice_score) <- rownames(datatoannotdup@meta.data)
    
    matrixscoreforHeatmap <-  SeuratObject::CreateSeuratObject(counts = t(matrice_score),meta.data = datatoannotdup@meta.data)
    
    matrixscoreforHeatmap$RNA@scale.data <- t(matrice_score)
    
    print(
    Seurat::DoHeatmap(object = matrixscoreforHeatmap,features = colnames(matrice_score) ,group.by = "seurat_clusters",size = 7)+
    ggplot2::theme(text = element_text(size = 20))
    )
    
    #suppression de datatoannotdup
    rm(datatoannotdup)
    
} # accolade fermante de la fonction cluster_annotated_SC_main
