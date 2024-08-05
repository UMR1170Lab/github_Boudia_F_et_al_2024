################################LibrarySizeBarplot###################################
# packages requis

# Cette fonction nécessite de charger les packages purrr et  rlang.

########### liste des arguments ############ 
# Data, dataframe avec en colonne les individus et en ligne les variables.
# Metadata : dataframe avec les métadatas (colonnes de métadatas). Le nombre de lignes de Metadata doit être égal au nombre de colonnes de Data. 
# ColMetadata : colonne du tableau de métadatas. Ne pas mettre entre guillemets le nom de la colonne.
#Couleur : vecteur de couleur (spécifiée entre guillemets comme des chaînes de caractères). Ce vecteur doit avoir une longueur égale au nombre de niveaux de la colonne de Metadata sélectionnée (argument ColMetadata) considérée comme un facteur. Pour choisir la couleur associée à chaque niveau de ColMetadata, il faut spécifié les couleurs dans l'ordre où les niveaux du facteur associé à ColMetadata apparaissent quand on exécute la commande du genre levels(factor(Colmetadata),exclude = NULL)
# xlabtitle : chaîne de caractères entre guillemets avec le nom de l'axe des x.
#ylabtitle : chaîne de caractères entre guillemets avec le nom de l'axe des y.
# TitreFigure : chaîne de caractères entre guillemets avec le titre de la figure. 
# TransLog (logical) : soit TRUE ou FALSE. Si TRUE, les données sont transformées (log(Data + 1)) sinon elles ne sont pas transformées pour la représentation.
LibrarySizeBarplot <- function(Data,Metadata,ColMetadata,Couleur,xlabtitle,ylabtitle,TitreFigure,TransLog){
  
  #traitement de l'argument ColMetadata
  
  ## traitement avec rlang
  
  ColMetadata_gui <- enexpr(ColMetadata)# l'objet créé va être traité comme une chaîne de caractères avec guillemets.
  
  
  # création d'un tableau qui associe à chaque niveau du facteur associé à la colonne ColMetadata
  
  AssoColorCondition <- data.frame(couleur = Couleur,cond = levels(factor(Metadata[,which(colnames(Metadata) == ColMetadata_gui)],exclude = NULL)))
  
  
  
  # création du vecteur de couleurs qui va être utilisé pour le barplot
  
  ColorBarplot <- sapply(as.character(Metadata[,which(colnames(Metadata)==ColMetadata_gui)]), function(x){
    
    return(AssoColorCondition$couleur[which(AssoColorCondition$cond == x)])
    
  })
  
  print(ColorBarplot)
  print(class(ColorBarplot))
  # représentation en barplot de la taille des librairies
  if(TransLog == FALSE){
    
    ## cas où les données ne sont pas loggées.
    librarysizebarplot <- barplot(unname(colSums(Data)), xlab = paste0(xlabtitle), ylab = paste0(ylabtitle), main = paste0(TitreFigure),col = as.character(ColorBarplot),border = as.character(ColorBarplot), xaxt = "n")
    
  }else{
    
    ## cas où les données sont loggés
    
    
    librarysizebarplot <- barplot(unname(colSums(log(Data+1))), xlab = paste0(xlabtitle), ylab = paste0(ylabtitle), main = paste0(TitreFigure),col = as.character(ColorBarplot),border = as.character(ColorBarplot), xaxt = "n")
    
  }
  
}
