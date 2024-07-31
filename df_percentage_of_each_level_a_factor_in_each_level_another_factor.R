######################## fonction qui donne la composition en niveaux d'un facteur dans chacuns des niveaux d'un autre facteur ######################

#' @title df_percentage_of_each_level_a_factor_in_each_level_another_factor



df_percentage_of_each_level_a_factor_in_each_level_another_factor <- function(df,statcol,otherfaccol){
    
    #On génère toutes les combinaisons possibles des niveaux des 2 facteurs d'intérêt
    comb <- base::expand.grid( unique(df[[statcol]]),unique(df[[otherfaccol]]))
    
    colnames(comb) <- c(statcol,otherfaccol)
    
    #On compte le nombre de lignes du data.frame df qui correspondent à chacune des combinaisons possibles des deux facteurs d'intérêts, on stocke le résultat dans un vecteur qui correspond à une colonne nommée effectif du dataframe effectifcombi qui correspond au dataframe comb auquel on a ajouté la colonne effectif.
    
    effectifcombi <- data.frame(
    comb,
    effectif = unlist(lapply(1:nrow(comb),function(z){
        
        print(df[which(df[[otherfaccol]]==comb[[otherfaccol]][z] & comb[[statcol]] == df[[statcol]][z]),])
        
        return(nrow(df[which(df[[otherfaccol]]==comb[[otherfaccol]][z] & df[[statcol]] == comb[[statcol]][z]),]))
        
    }))
    )
    
    
    effectifcombi <- data.frame(effectifcombi,
    percentage = unlist(lapply(1:nrow(effectifcombi),function(z){
        
        return(
        ((effectifcombi$effectif[z])/(sum(effectifcombi$effectif[which(effectifcombi[[statcol]]==effectifcombi[[statcol]][z])])))*100
        )
        
        
        
    }))
    
    )
    
    return(effectifcombi)
    
    
    
}
