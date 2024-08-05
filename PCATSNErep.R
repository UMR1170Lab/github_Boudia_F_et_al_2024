######################################## PCATSNErep ######################################

## fonction qui permet d'avoir directement les représentations en PCA et en TSNE.

## Packages

# Cette fonction nécessite de charger les package purrr, rlang et tidyverse. Elle nécessite également le chargement des packages FactoMineR et factoextra ainsi que Rtsne. Nécessite aussi magrittr

## Arguments

# Data, un objet de classe matrice ou dataframe avec en colonne les variables et en lignes les samples (individus, échantillons, cellules ...)

# Metadata : un tableau de Metadata avec en colonnes les différents type de métadonnées (ex sexe, phénotype ...) pour chaque individus de Data (le nombre de lignes de Metadata doit être égal au nombre de lignes de Data). Les lignes associées à chaque individu doivent être rangées dans le me ordre dans Data et dans Metadata.

# ColMetadata : la colonne de Metadata qui va être le facteur déterminant dans quelle couleur va être le point associé à chaque individu (NOM A NE PAS METTRE ENTRE GUILLEMETS).

# ColorPoint : vecteur de couleurs spécifiant pour chaque niveau du facteur correspondant à la colonne ColMetadata de Metadata à quelle couleur celui-ci est associé (détermine la couleur des points selon l'appartenance de chacuns des individus à un des niveaux de ColMetadata). L'ordre des couleurs est associé à l'ordre dans lequel les niveaux du facteur ColMetadata

#PCA_Arg : un objet de classe list qui va contenir les arguments choisis pour faire la PCA. Les différents éléments de la liste qui servent pour faire la PCA sont :
# scaling,ndim, DimX, DimY,titlePCA. Les noms des éléments doivent être EXPLICITEMENT DONNE sinon la fonction ne va pas marcher.
# scaling,ndim,DimX = Dim.1,DimY = Dim.2,title
# exemple correct : PCA_Arg = list(scaling = TRUE,ndim = 2, DimX = "Dim.1",DimY = "Dim.2",title = "PCA")
# L'ordre dans lequel sont donné les éléments de la liste n'a pas d'importance.
# Les dimensions doivent être spécifiées sous la forme "Dim.N" (chaîne de caractères entre guillemets) où N est le nombre entier spécifiant la dimension
# Les dimensions de la PCA que l'on souhaite représenter doivent être inférieures ou égales à ndim
# ndim doit être supérieure ou égale à max()
# dims = 2, perplexity = 20, max_iter = 6000
#TSNE_arg : un objet de class list qui va contenir les éléments suivants :
# ndim : le nombre de dimensions calculées par la fonction Rtsne
# perp : la perplexité.
#nb_itermax : le nombre d'itérations maximales.



PCATSNErep <- function(Data,Metadata,ColMetadata,PCA_Arg,TSNE_arg,ColorPoint){
  
  # traitement de la variable ColMetadata avec rlang.
  ColMetadata_gg <- rlang::enquo(ColMetadata)
  
  ColMetadata_gui <- enexpr(ColMetadata)# l'objet créé va être traité comme une chaîne de caractères avec guillemets.
  ColMetadata_Char <- colnames(Metadata)[which(colnames(Metadata)==ColMetadata_gui)]
  
  # spécification des dimensions de la PCA à représenter en abscisse et en ordonnées
  # DimAbsPCA <- PCA_arg[["DimX"]]
  # 
  # DimOrdPCA <- PCA_arg[["DimY"]]
  
  
  data_to_represent <- invoke_map(list(function(x,scaling = PCA_Arg[["scaling"]],ndim = PCA_Arg[["ndim"]],DimX = PCA_Arg[["DimX"]],DimY = PCA_Arg[["DimY"]],colmetadata = ColMetadata_Char,metadata= Metadata){
    ##########################début de la fonction qui prépare les données pour la représentation de la pca.##############
    
    sourcePCA <- PCA(x,scale.unit = scaling,ncp = ndim,graph =FALSE)#réalisation de la PCA
    coord_ind_activ_pca <- get_pca_ind(sourcePCA)
    coord_ind_activ_pca <- coord_ind_activ_pca$coord %>% as.data.frame() %>% mutate(category = metadata[,which(colnames(metadata)==colmetadata)]) %>% as.data.frame()#récupération des coordonnées des points de la PCA et ajout de la colonne de metadata.
    
    
    colnames(coord_ind_activ_pca)[which(colnames(coord_ind_activ_pca)=="category")] <- colmetadata#on renomme la colonne de metadata avec le vrai nom de la variable de metadata.
    
    
    EigenValuePCA <- as.data.frame(get_eigenvalue(sourcePCA))
    # pourcentage de variance expliquée par l'axe des abscisse
    PerVar1 <- EigenValuePCA[which(rownames(EigenValuePCA)==DimX),2]
    # print(paste0(rows1," : ",PerVar1," %"))
    #Pourcentage de variance expliquée par l'axe des ordonnées
    PerVar2 <- EigenValuePCA[which(rownames(EigenValuePCA)==DimY),2]
    
    return(list(DataToRep = coord_ind_activ_pca, TextLabel = c(paste0(DimX, " = ",PerVar1, " %"),paste0(DimY, " = ",PerVar2," %")),selecCol = c(DimX,DimY)))
    
    
    ############ fin de la fonction qui prépare les données pour la représentation de la pca.  
  },function(x,ndim = TSNE_arg[["ndim"]],perp = TSNE_arg[["perp"]],nb_itermax = TSNE_arg[["nb_itermax"]],colmetadata = ColMetadata_Char,metadata= Metadata){
    
    ###################### début de la fonction qui prépare les données pour la représentation de la tsne ######################   
    
    tsne_list <- Rtsne(x,dims = ndim, perplexity = perp, max_iter = nb_itermax)# réalisation de la TSNE
    
    
    
    tsne_coord <- as.data.frame(tsne_list$Y) %>% mutate(category = metadata[,which(colnames(metadata)==colmetadata)])#récupération des coordonnées des points de la PCA et ajout de la colonne de metadata.# obtention d'un dataframe avec les coordonnées issues de la TSNE pour chaque individu
    
    colnames(tsne_coord)[which(colnames(tsne_coord)=="category")] <- colmetadata#on renomme la colonne de metadata avec le vrai nom de la variable de metadata.
    
    #################### probablement qu'il faut que je nomme les éléments de la liste #################################
    return(list(DataToRep = tsne_coord, TextLabel = c("tsne-1","tsne-2"),selecCol = c("V1","V2")))
    
    #fin de la fonction qui génère les données pour la représentation de la tsne
    
  }
  )############# fin de la liste de fonctions
  ,x = Data
  )################################fin de la première fonction invoke_map ######################################
  
  # return(data_to_represent)
  
  
  invoke_map(list(function(DataToRep,TextLabel,selecCol){
    
    
    
    print(ggplot(DataToRep,aes_string(x = selecCol[1],y= selecCol[2],colour = ColMetadata_gg))+
            geom_point()+
            scale_color_manual(values = ColorPoint)+
            xlab(TextLabel[1])+
            ylab(TextLabel[2])+
            theme_bw()
    )
    
  },function(DataToRep,TextLabel,selecCol){
    
    print(ggplot(DataToRep,aes_string(x = selecCol[1],y = selecCol[2],colour = ColMetadata_gg))+
            geom_point()+
            scale_color_manual(values = ColorPoint)+
            xlab(TextLabel[1])+
            ylab(TextLabel[2])+
            theme_bw()
          
          
    )
    
    
  }
  
  
  
  
  
  ),data_to_represent)
  #   
  
}
