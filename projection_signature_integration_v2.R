#Objectif de la fonction : projeter les statuts de signature génique sur les intégrations et voir systématiquement. Représentations barplot comportement des signatures au sein de chaque cluster.
#integrated_dataset : objet Seurat contenant une intégration de différents batchs de Single-Cell.
#norm_batch_dataset : une list R. Chaque élément de cette list R correspond à un objet Seurat qui contient la matrice normalisée d'un des batchs qui est inclu dans #l'intégration. Pour chaque batch contenu dans l'intégration,un unique élément correspondant à l'objet Seurat normalisé qui lui est associé doit obligatoirement être #présent dans cette list. Les noms des éléments de cette liste (indication du batch, du dataset d'origine des cellules) doivent correspondre aux niveaux des facteurs dans #la colonne orig.ident du data.frame meta.data de l'objet Seurat de l'intégration.
#interest_biological_signature : vecteur de class character ayant comme éléments les gene symbol de l'ensemble des gènes qui constituent une signature génique.
#name_biological_signature : un élément de class character qui correspond au nom que l'on veut donner à la signature biologique.
#path_signature_statut_cells_save_rds (élément de chaîne caractère) path et nom du fichier rds qui va être généré une fois que le statut des cellules de l'intégration #concernant la signature d'intérêt aura été calculé. Si ce fichier a déjà été généré auparavant lors d'un précédent usage de la fonction, alors c'est ces résultats qui #sont utilisés et la fonction  IndiceRepSignaturePval_v5 n'est pas utilisé de nouveau (gain de temps qui peut être important).
#minimal_prop_cells : (entre 0 et 1 : par exemple si on veut 2% comme seuil, on met 0.02) proportion minimale de cellules de la matrice d'expression dans laquelle doit être détecté un gène pour qu'il soit conservé pour la suite des analyses.
#'number_random_signatures (integer) défaut de 100. Nombre de signatures aléatoires pour lequel on calcule le taux de détection (pourcentage de gènes détectés)
#
#
#threshold_pval (numeric de une valeur) : threshold de p-value pour considérer une cellule comme possédant plus de gènes de l'ensemble de gènes testés que ce que l'on attendrait par hasard.
#used_random_seed (integer de une valeur) : graine aléatoire pour garantir la reproductibilité des résultats
#path_pdf_figures (élément de class character) : path et nom que l'on donne au pdf qui va contenir les figures générées par la fonction.
# width_pdf : (élément de class numeric) : largeur du pdf
# height_pdf : (élément de class numeric) : hauteur du pdf.
projection_signature_integration_v2 <- function(integrated_dataset,norm_batch_dataset,interest_biological_signature,path_signature_statut_cells_save_rds,name_biological_signature,minimal_prop_cells,number_random_signa,threshold_pval,used_random_seed,path_pdf_figures,width_pdf,height_pdf){#accolade ouvrante de la fonction projection_signature_integration_v1
    
    #chargement des packages
    
    require(Seurat)
    require(readr)
    require(purrr)
    
    
    # 1) Dans un premier temps, on calcule le statut de chaque cellule dans les datasets séparés.
    
    print(interest_biological_signature)
    print(file.exists(path_signature_statut_cells_save_rds))
    if(file.exists(path_signature_statut_cells_save_rds)==FALSE){ #accolade ouvrante : instructions si le fichier indiqué au path path_signature_statut_cells_save_rds #n'existe pas déjà.
        print("Calcul de signatures pour chacuns des dataframes individuels")
        cell_signature_batch_statut <- purrr::pmap(.l= list(data = norm_batch_dataset, nom = as.list(names(norm_batch_dataset))),.f = function(data,nom){
            #accolade ouvrante de la fonction du pmap
            print(nom)
            # #calcul du statut des cellules dans chacun des jeux de données (calcul qui est réalisé de façon individuel)
            IndiceRepSignaturePval_v5(data@assays$RNA@counts,
            biological_signature = interest_biological_signature,
            min_prop_cells_detected = minimal_prop_cells,
            number_random_signatures = number_random_signa,
            thresholdpval = threshold_pval,
            random_seed = used_random_seed)
            
            # #accolade fermante de la fonction pmap.
        })
        
        print("exportation des résultats pour chacun des jeux de données")
        #exportations des résultats de la fonction IndiceRepSignaturePval_v5
        readr::write_rds(cell_signature_batch_statut,path_signature_statut_cells_save_rds)
        #accolade fermante : instructions si le fichier indiqué au path path_signature_statut_cells_save_rds n'existe pas déjà.
    }
    else{#accolade ouvrante : instructions si le fichier indiqué au path path_signature_statut_cells_save_rds existe.
        
        cell_signature_batch_statut <- readr::read_rds(path_signature_statut_cells_save_rds)
        #accolade fermante  : instructions si le fichier indiqué au path path_signature_statut_cells_save_rds existe.
    }
    #2) On ne récupère que le statut des cellules au regard de la signature biologique d'intérêt dans la list de résultats générés plus haut.
    cells_statut_pvalue <- lapply(cell_signature_batch_statut,function(z){
        return(
        
        z[["statut_pvalue"]]
        
        )
    })
    #3)On génère le dataframe qui permet de savoir à quel chiffre est associé chaque batch dans le dataframe de métadatas de l'intégration.
    association <- association_number_condition(integrated_dataset)
    #4)Dans la list qui donne pour chaque batch de l'intégration, le statut des cellules concernant la signature biologique d'intérêt, on ajoute à l'identifiant de chaque #cellule un chiffre précédé d'un _ qui indique à quel batch appartient la cellule. Pour cela, on utilise le dataframe association généré ci-dessus.
    for(i in 1:length(cells_statut_pvalue)){#accolade ouvrante de la boucle
        
        names(cells_statut_pvalue[[i]]) <- paste0(names(cells_statut_pvalue[[i]]),"_",association$number[which(association$condition == names(cells_statut_pvalue)[i])])
        
        
        
        #ci-dessous, accolade fermante de la boucle.
    }
    #5) On génère un vecteur en fusionnant l'ensemble des éléments de la list cells_statut_pvalue
    cells_statut_pvalue_concat <- cells_statut_pvalue %>% do.call(c,.)
    #6) On vire la première partie des identifiants des cellules (noms des éléments du vecteur) qui se trouve avant le premier . dans les identifiants. Cette partie #correspond au noms du batch dont proviennent chacune des cellules.
    names(cells_statut_pvalue_concat) <- names(cells_statut_pvalue_concat) %>%  strsplit("[.]") %>% lapply(function(x){return(paste0(x[-1],collapse="."))})
    #7) On génère une colonne dans le dataframe de metadatas de l'intégration qui donne pour chaque cellule son statut concernant la signature biologique d'intérêt.
    integrated_dataset@meta.data[[name_biological_signature]] <- (lapply(rownames(integrated_dataset@meta.data),function(z){
        
        return(cells_statut_pvalue_concat[which(names(cells_statut_pvalue_concat)==z)])
        
        
    }) %>% unlist())
    # Les étapes suivantes permettent de générer des représentations graphiques.
    pdf(path_pdf_figures,width = width_pdf,height = height_pdf)
    #8) Représentation Seurat. Toutes les cellules (statut) Les cellules de tous les batchs sont représentées.
    print(
    Seurat::DimPlot(integrated_dataset,
    reduction="umap",
    group.by = name_biological_signature,
    cols= c("NOT_SIGNI"="grey",
    "POS_SIGNI" = "red",
    "NEG_SIGNI" = "blue"),
    order = c("POS_SIGNI","NEG_SIGNI","NOT_SIGNI"))+
    ggplot2::labs(title=paste0(name_biological_signature," All cells of the integration"),subtitle = paste0("min_prop_cells :", minimal_prop_cells," number_random_signa : ",number_random_signa," thresh pval :",threshold_pval," random seed used : ",used_random_seed))
    )
    
    # représentation du pourcentage de pos signi, neg signi et non signi pour la signature d intérêt dans chacun des clusters.
    df_rep_of_each_level_a_factor_in_each_level_another_factor(integrated_dataset@meta.data,"seurat_clusters",name_biological_signature,"percentage",personnal_colors = TRUE,colchoose=c("NOT_SIGNI"="grey","POS_SIGNI" = "red","NEG_SIGNI" = "blue"))
    
    #9) Pareil qu'en 8 sauf que  cette fois-ci, dans la umap, on ne présenrve que les points qui appartiennent à un seul batch et on ne regarde le statut que pour les #cellules de ce batch. On réitère l'opération pour toutes les cellules du batch.
    lapply(unique(integrated_dataset@meta.data$orig.ident),function(z){#accolade ouvrante de la fonction du lapply.
        
        
        #On crée une colonne spécifique pour chaque batch.
        # Dans ces colonnes, les indications concernant le statut des cellules pour la signature d'intérêt ne figure que pour le batch concerné. On remplace les indications
        #concernant les cellules des autres batchs par others_condition.
        
        integrated_dataset@meta.data[[paste0(name_biological_signature,"_",z)]] <- integrated_dataset@meta.data[[name_biological_signature]][which(integrated_dataset@meta.data$orig.ident != z)] <- "others_condition"
        
        #représentation Seurat avec les indications ne figurant que pour un seul batch.
        print(
        Seurat::DimPlot(integrated_dataset,
        reduction="umap",
        group.by = name_biological_signature,
        cols= c("NOT_SIGNI"="grey",
        "POS_SIGNI" = "red",
        "NEG_SIGNI" = "blue",
        "others_condition"="white"),
        order = c("POS_SIGNI","NEG_SIGNI","NOT_SIGNI","others_condition"))+
        ggplot2::labs(title=paste0(name_biological_signature," ",z),subtitle = paste0("min_prop_cells :", minimal_prop_cells," number_random_signa : ",number_random_signa," thresh pval :",threshold_pval," random seed used : ",used_random_seed))
        
    
        
        
        
        )
        

        
        
        
        #
        
        #ci-dessous, accolade fermante de la fonction du lapply  
    }
    )
    dev.off()
    #accolade fermante de la fonction projection_signature_integration_v2
}
