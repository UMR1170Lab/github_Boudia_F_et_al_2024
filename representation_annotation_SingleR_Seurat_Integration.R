@title representation_annotation_SingleR_Seurat_Integration




representation_annotation_SingleR_Seurat_Integration <- function(normalized_data_list,Seurat_Integration,reference_SingleR, vector_col_annotation = NULL,pathpdf,height_pdf,width_pdf,pt_size_dimplot,path_rds_file){
    # ci-dessus accolade ouvrante de la fonction representation_annotation_SingleR_Seurat_Integration
    
    require(readr)
    require(SingleR)
    
    # On annote chacun des datasets dont provient l intégration séparément.
    annot_all_normalized_data <- lapply(normalized_data_list,function(z){#accolade ouvrante de la fonction du lapply
        print("Annotation SingleR")
        #Utilisation de la fonction annotation_SingleR_dataframe_output
        annotation_SingleR_dataframe_output(z,reference = reference_SingleR)
        
        #accolade fermante de la fonction du lapply
    })
    
    #On renomme les éléments de la list contenant les annotations.
    print("On renomme les éléments de la list générée")
    names(annot_all_normalized_data) <- names(normalized_data_list)
    
    # Seurat_Integration2 correspond à l intégration Seurat.
    Seurat_Integration2 <- Seurat_Integration
    
    #On génère le dataframe qui permet de savoir à quel chiffre est associé chaque batch dans le dataframe de métadatas de l'intégration.
    association <- association_number_condition(Seurat_Integration2)
    
    
    #Boucle pour renommer les lignes des df d'annotations des cellules de chacuns des datasets qui a été employé pour l intégration.
    for(i in 1:length(annot_all_normalized_data)){#accolade ouvrante de la boucle
        #renommage des lignes en ajoutant le numéro du dataset dans l intégration précédé d'un _ au nom de la ligne préexistant.
        rownames(annot_all_normalized_data[[i]]) <- paste0(rownames(annot_all_normalized_data[[i]]),"_",association$number[which(association$condition == names(annot_all_normalized_data)[i])])
        #ci-dessous, accolade fermante de la boucle.
    }
    
    #récupération des sous-df qui contiennent uniquement la colonne de données qui  correspond au type de cellule inféré par SingleR en
    # utilisant le dataset de Novershtern.
    # Fusion des différents dfs ainsi générés
    ann_all <- lapply(annot_all_normalized_data,function(z){return(z[,which(colnames(z)=="all_single_R_annotation"),FALSE])}) %>% {do.call(rbind,.)}
    
    # On transforme le noms des lignes du df précédemment généré pour qu'elles correspondent exactement aux noms des lignes du df de meta.data # de l intégration Seurat.
    
    rownames(ann_all) <-  (strsplit(rownames(ann_all),".",fixed=TRUE) %>% lapply(function(z){return(z[2])}) %>% unlist())
    
    #Génération d un df qui correspond au df de meta.data avec une colonne d annotation SingleR en plus.
    integration_data_annotated <- base::merge(Seurat_Integration2@meta.data,ann_all,"row.names")
    
    
    
    #Le nouveau df de meta.data de l intégration correspond au df avec les annotations.
    Seurat_Integration2@meta.data <- integration_data_annotated
    
    # On remplace les NAs pour les cellules pour lesquelles SingleR n'a pas pas réussi à attribuer un type cellulaire par undetermined
    Seurat_Integration2@meta.data$all_single_R_annotation[which(is.na(Seurat_Integration2@meta.data$all_single_R_annotation)==TRUE)] <- "undetermined"
    
    
    #les noms des cellules vont correspondre aux noms des cellules dans le df Seurat_Integration2
    rownames(Seurat_Integration2@meta.data) <- Seurat_Integration2@meta.data$Row.names
    
    print("premières lignes de Seurat_Integration2 metadata")
    
    print(class(Seurat_Integration2))
    
    print(head(Seurat_Integration2@meta.data))
    
    #On représente en umap les cellules colorés en fct du type cellulaire attribué par SingleR.
    
    if(is.null(vector_col_annotation)==FALSE){
        #ci-dessus accolade ouvrante instructions si vecteur de couleurs défini
        
        print("vecteur de couleurs d annotation fourni")
        
        print("début de la génération de la représentation")
        
        pdf(pathpdf,height = height_pdf,width = width_pdf)
        
        print(
        Seurat::DimPlot(Seurat_Integration2,reduction="umap",group.by = "all_single_R_annotation",cols = vector_col_annotation,pt.size = pt_size_dimplot)
        )
        
        dev.off()
        print("fin de la génération de la représentation")
        #ci-dessus accolade ouvrante instructions si vecteur de couleurs défini
    }else{
        #ci-dessus accolade ouvrante instructions si vecteur de couleurs non défini
        
        print("vecteur de couleurs d annotation non fournie")
        print("Génération de la représentation")
        pdf(pathpdf,height = height_pdf,width = width_pdf)
        
        print(
        Seurat::DimPlot(Seurat_Integration2,reduction = "umap",group.by = "all_single_R_annotation", pt.size = pt_size_dimplot)
        )
        
        dev.off()
        
        print("fin de la génération de la représentation")
        
        #ci-dessous accolade ouvrante instructions si vecteur de couleurs non défini
    }
    
    # exportation au format rds du fichier de metadata contenant les annotations SingleR Novershtern
    print("exportation du dataframe contenant les annotations SingleR")
    print(path_rds_file)
    readr::write_rds(Seurat_Integration2@meta.data,path_rds_file)
    
    # ci-dessous accolade fermante de la fonction representation_annotation_SingleR_Seurat_Integration
}
