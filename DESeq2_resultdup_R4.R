#analyses différentielles pour ce RNAseq
DEseq2_resultdup_R4 <- function(txiFile,metadataTab,threshold_keep_gene=10){
  
  # fonction pmap pour répéter la même opération pour chaque paire de conditions comparés 
  resultExpDiff <- purrr::pmap(.l = list(
    x1 = as.list(combn(unique(metadataTab$condition),2)[1,]),# premier élément des paires de #conditions à comparer.
    x2 = as.list(combn(unique(metadataTab$condition),2)[2,])),# second élément des paires de#conditions à comparer.
    .f = function(x1,x2){# fonction à appliquer
      
      ##############création d'un objet DESeqDataSet nommé dds #########
      dds <- DESeqDataSetFromTximport(
        txiFile,
        metadataTab,# sélection de la partie du dataframe de métadatas associée aux conditions comparées 
        ~condition)
      
      #############Filtrage de l'objet dds #################
    
      #On ne garde que les gènes qui ont au moins 10 counts pour l'analyse.
      keep <- rowSums(counts(dds)) >= threshold_keep_gene
      dds <- dds[keep,]
      
      rm(keep)#suppression de l'objet keep une fois utilisé pour sélectionner les lignes à conserver : voir ci-dessus
      
      ################ fixation du facteur de référence #################################
      
      #les facteurs de référence sont les conditions figurant dans la première list de la list donnée en entrée de dds.
      dds$condition <- relevel(dds$condition,ref = paste0(x1))
      
      ############## Réalisation de l'analyse différentielle#####################################
      
      dds <- DESeq(dds)# fonction qui effectue l'analyse différentielle
      
      #################  résultats de l'analyse différentielle qui vont être renvoyés par la fonction de pmap #####################
      
      result_paire <- list(sansShrinkage = data.frame(DESeq2::results(dds,name= paste0("condition_",x2,"_vs_",x1)),gene_id = rownames(results(dds,name= paste0("condition_",x2,"_vs_",x1))))  %>% arrange(padj) %>% conv_gene_tab(genesconvTab = genes_conv),
                           avecShrinkage = data.frame(lfcShrink(dds,coef=paste0("condition_",x2,"_vs_",x1),type="apeglm"),
                                                       gene_id = rownames(lfcShrink(dds,coef=paste0("condition_",x2,"_vs_",x1),type="apeglm"))) %>% dplyr::arrange(padj) %>% conv_gene_tab(genesconvTab = genes_conv)) 
      # on donne les résultats avec et sans le shrinkage.
      
      
      ############################ suppression du dds ###########################################
      
      rm(dds)
      
      return(result_paire)
      
      
    }# accolade fin de la fonction de pmap
  )# parenthèse fin de la fonction de pmap
  
  
  
  ###########################on renvoie le résultat de la fonction, la liste resultExpDiff ######################################
  names(resultExpDiff) <- paste0(combn(unique(metadataTab$condition),2)[2,],"_vs_",combn(unique(metadataTab$condition),2)[1,])
  
  
  
  return(resultExpDiff)  
  
  
  
}

