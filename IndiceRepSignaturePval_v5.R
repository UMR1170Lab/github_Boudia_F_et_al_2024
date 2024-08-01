#Cette fonction peut être utiliser plusieurs fois en faisant varier la graine aléatoire de façon à pouvoir calculer des moyennes des p-values
IndiceRepSignaturePval_v5 <- function(exprMat,biological_signature,min_prop_cells_detected,number_random_signatures = 100,thresholdpval,random_seed){
#accolade ouvrante fonction IndiceRepSignaturePval_v5
  
#Fixation d'une graine aléatoire afin d'assurer la reproductibilité des données.  
print("Fixation d'une graine aléatoire afin d'assurer la reproductibilité des données.")
set.seed(random_seed) 
 
print("Sélection des gènes qui vont être utilisé pour tester les signatures. La matrice résultante ne contenant que ces gènes sera nommé used_exprMat")

#Sélection des gènes qui vont être utilisé pour tester les signatures. La matrice résultante ne contenant que ces gènes sera nommé used_exprMat
used_exprMat <- exprMat[which((apply(exprMat,1,function(x){length(which(x>0))}) %>% unlist())>= min_prop_cells_detected*ncol(exprMat)),]
######################################### détermination des taux de détection pour les signatures aléatoires #################################### 
# Génération de 100 signatures aléatoires (en fait on garde que la position des gènes dans la matrice)
echantillonnage <- replicate(number_random_signatures,sample(rownames(used_exprMat),length(intersect(rownames(used_exprMat),biological_signature))),simplify = FALSE) %>% lapply(function(x){
  #accolade ouvrante de la fonction permettant de récupérer la position des gènes
  
  which(rownames(used_exprMat) %in% x)
  
  
  
#accolade fermante de la fonction  
})

print("pour chacune des goutelettes, on applique pour chacune des signatures.")
test <- apply(used_exprMat,2,function(z){#pour chacune des goutelettes, on applique pour chacune des signatures.
  
  
droplets <- z
  
 percentages <- lapply(echantillonnage,function(y){#établissement du pourcentage de nombre de gènes détectés pour chacune des signatures aléatoires.
  
   (length(which(droplets[y]>0))/length(droplets[y]))*100
   
   
   
 }# établissement du pourcentage pour chacune des signatures aléatoires. 
 ) 
  
  
  return(percentages)
  
  
})

print("pour chacune des droplets, on stocke dans un vecteur les pourcentages de gènes détectés pour chacune des signatures aléatoires")
#pour chacune des droplets, on stocke dans un vecteur les pourcentages de gènes détectés pour chacune des signatures aléatoires
test <- lapply(test,unlist)
#les noms des éléments de la liste de stockage vont correspondre aux noms de la colonne de la matrice.

print("Les noms des éléments de la liste de stockage vont correspondre aux noms de la colonne de la matrice.")
names(test) <- colnames(used_exprMat)
#pour chacune des cellules, on détermine la moyenne du taux de détection des signatures aléatoires (en pourcentage)
print("pour chacune des cellules, on détermine la moyenne du taux de détection des signatures aléatoires (en pourcentage)")
mean_random_detection_rate_per_cell <- lapply(test,mean) %>% unlist()
##################### Détermination pour la véritable signature biologique des taux de détection associés à chacune des droplets. #############

print("Détermination pour la véritable signature biologique des taux de détection associés à chacune des droplets.")
detection_rate_signature <- apply(used_exprMat,2,function(z){
  
    
  return(
  (length(which(z[which(rownames(used_exprMat) %in% intersect(rownames(used_exprMat),biological_signature))]>0))/length(z[which(rownames(used_exprMat) %in% intersect(rownames(used_exprMat),biological_signature))]))*100
  )
})
#Calcul du ratio detection_rate_signature sur mean_random_detection_rate_per_cell 
print("Calcul du ratio detection_rate_signature sur mean_random_detection_rate_per_cell") 
ratio_interest_sign_vs_random <- (detection_rate_signature/mean_random_detection_rate_per_cell)
### calcul des p-values : probabilité d'obtenir un taux de détection donné sous l'hypothèse que la signature a été générée aléatoirement.

print("calcul des p-values : probabilité d'obtenir un taux de détection donné sous l'hypothèse que la signature a été générée aléatoirement.")
pvaltoret <- purrr::pmap(.l = list(random = test , signature = detection_rate_signature),.f = function(random,signature){
  
return(pvaldistri_v3(random,signature))
  
}) %>% unlist()
####################################### Détermination du statut d'enrichissement de la signature pour chacune des cellules #######################

print("Détermination du statut d'enrichissement de la signature pour chacune des cellules")
statutpval <- dplyr::case_when(pvaltoret < 0 & abs(pvaltoret) < thresholdpval ~ "NEG_SIGNI",
                                                                pvaltoret < 0 & abs(pvaltoret) >= thresholdpval ~ "NOT_SIGNI",
                                                                pvaltoret > 0 & abs(pvaltoret) < thresholdpval ~ "POS_SIGNI",
                                                                pvaltoret > 0 & abs(pvaltoret) >= thresholdpval ~ "NOT_SIGNI")
#On nomme les éléments du vecteur qui donne le statut associé à chaque goutelette en fonction de sa p-value
names(statutpval) <- names(pvaltoret)
######################################### On renvoie une liste contenant pour chacune des droplets les valeurs des paramètres qui permettent de juger de l'enrichissement d'une signature dans une cellule donnée ## 
    return(list("ratio_interest_sign_vs_random" = ratio_interest_sign_vs_random,
               "detection_rate_signature" = detection_rate_signature,
               "mean_random_detection_rate_per_cell" = mean_random_detection_rate_per_cell,
               "p-value" = pvaltoret,
               "abs_pvalue" = abs(pvaltoret), 
               "statut_pvalue" = statutpval            
              ))
  
##accolade fermante de la fonction IndiceRepSignaturePval_v5
}

