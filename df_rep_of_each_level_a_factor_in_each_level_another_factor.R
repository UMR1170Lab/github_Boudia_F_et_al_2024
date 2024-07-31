################## df_rep_percentage_of_each_level_a_factor_in_each_level_another_factor


######################### fonction qui va être utiliser pour faire les représentations en ggplot de la composition des différents clusters


#' @title df_rep_of_each_level_a_factor_in_each_level_another_factor

df_rep_of_each_level_a_factor_in_each_level_another_factor <- function(df,statcol,otherfaccol,optionrep,personnal_colors=FALSE,colchoose = c("blue","red")){
    
    if(class(optionrep)!="character"){
        
        stop("The value given for the argument optionrep must be of class character (must be either brut or percentage)")
        
        
    }
    
    if(optionrep %in% c("brut","percentage")==FALSE){
        
        stop("The value given for the argument optionrep must be either the character chain brut or percentage (must be of class character)")
        
        
    }
    
    if(class(personnal_colors)!="logical"){
        
        stop("The value given for the argument personnal_colors must be of class logical")
        
    }
    
    
    require(ggplot2)
    require(chameleon)
    
    #On génère toutes les combinaisons possibles des niveaux des 2 facteurs d'intérêt
    comb <- base::expand.grid( unique(df[[statcol]]),unique(df[[otherfaccol]]))
    
    colnames(comb) <- c(statcol,otherfaccol)
    
    #On compte le nombre de lignes du data.frame df qui correspondent à chacune des combinaisons possibles des deux facteurs d'intérêts, on stocke le résultat dans un vecteur qui correspond à une colonne nommée effectif du dataframe effectifcombi qui correspond au dataframe comb auquel on a ajouté la colonne effectif.
    
    effectifcombi <- data.frame(
    comb,
    effectif = unlist(lapply(1:nrow(comb),function(z){
        
        # print(df[which(df[[otherfaccol]]==comb[[otherfaccol]][z] & comb[[statcol]] == df[[statcol]][z]),])
        
        return(nrow(df[which(df[[otherfaccol]]==comb[[otherfaccol]][z] & df[[statcol]] == comb[[statcol]][z]),]))
        
    }))
    )
    
    ############# calcul des pourcentages des différents niveaux du facteur otherfaccol dans chaque niveau de statcol. ###############################
    effectifcombi <- data.frame(effectifcombi,
    percentage = unlist(lapply(1:nrow(effectifcombi),function(z){
        
        return(
        ((effectifcombi$effectif[z])/(sum(effectifcombi$effectif[which(effectifcombi[[statcol]]==effectifcombi[[statcol]][z])])))*100
        )
        
        
        
    }))
    
    )
    
    
    if(personnal_colors==FALSE){
        
        col_otherfaccol <- chameleon::distinct_colors(length(unique(df[[otherfaccol]])))[["name"]]
        
    }else{
        
        if(length(colchoose) != length(unique(df[[otherfaccol]]))){
            
            stop(paste0("The number of colors in the vector colchoose must be equal to the number of levels of the factor ",otherfaccol, " in the dataframe given as value of the argument df "))
        }
        else{
            
            col_otherfaccol <- colchoose
            
        }
        
        
        
    }
    
    
    
    ################################################## représentation ggplot ###################################################
    
    
    effectifcombi[[otherfaccol]] <- as.character(effectifcombi[[otherfaccol]])
    
    if(optionrep=="percentage"){
        
        
        print(
        ggplot(effectifcombi,aes_string(x=statcol, y="percentage", fill=otherfaccol))+ geom_bar(stat="identity")+
        theme_bw()+
        coord_flip()+
        scale_fill_discrete(breaks = sort(unique(df[[otherfaccol]])))+
        scale_fill_manual(values= col_otherfaccol)
        
        )
    }else{
        
        print(
        ggplot(effectifcombi,aes_string(x=statcol, y="effectif", fill= otherfaccol))+ geom_bar(stat="identity")+
        theme_bw()+
        coord_flip()+
        scale_fill_discrete(breaks = sort(unique(df[[otherfaccol]])))+
        scale_fill_manual(values= col_otherfaccol)
        
        )
        
        
    }
    
}

