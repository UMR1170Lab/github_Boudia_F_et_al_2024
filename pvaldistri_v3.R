
#' @title pvaldistri_v3
#' @description calcule la p-value d'un taux de détection pour une signature biologique sous l'hypothèse nulle que celle-ci est équivalente à une signature aléatoire.
#' @param distridetecrate numeric vector vector of distribution rate of the set of random signature that hav the same size that the signature of interest.
#' @param detectionratesignature detection rate of the signature of interest
#' @return return a numeric value which represent the probability that this detection rate can be observed with a random signature.  
#fonction pour calculer les p-values sur les distributions 
#distridetecrate :(numeric) : vecteur des taux de distribution sur les ensembles de gènes aléatoires
#detectionratesignature (numeric) : valeur du taux de détection pour la signature d'intérêt (pourcentage)

pvaldistri_v3 <- function(distridetecrate,detectionratesignature){
  
  
  verysmall <- 0
  verylarge <- 200000
  
#Génération d'un vecteur de densité
X = c(verysmall, density(distridetecrate)$x, verylarge)
#somme cumulée de la fonction de densité
SUM = cumsum(density(distridetecrate)$y)
# Vecteur normalisé sur 1
Y = c( 0, SUM/max(SUM) , 1)
# Valeur x de la distribution pour laquelle y = 0.5 (moitié  de l'AUC de la fonction de distribution 0.5) 
middistribution <- approxfun(Y[-c(1,length(Y))],X[-c(1,length(X))])(0.5)
if(detectionratesignature >= middistribution){
  
  Yrev <- rev(Y)
  
  pval <- approxfun(X,Yrev)(detectionratesignature)
  if(pval==0){
    
    pval <- 10^-10
  }
  return(pval)
  
}
else{
    
  
  pval <- approxfun(X,Y)(detectionratesignature)
    if(pval==0){
    
    pval <- 10^-10
    
  }
  return(-pval)
  
  
}
  
  
  
}
