






#This function allows to explore integration generated using the RPCA method provided by Seurat. The user can test a range of values both for the number of PCA dimensions to keep and for the resolution of the clustering and finally select the wanted values for these two parameters.



representation_summary_integration_Seurat_datasets_rpca_v1 <- function(rootdir_output_function,
random_seed,
vec_dimension,vec_resolution,
integrated_dataset_rds,color_vector_dataset,
gene_pop,normalized_data,
nom_reduction,equivalence_nom_sample_nom_figure,
ref_SingleR,vector_col_annota,
pt_size_dimplot_SingleR,
Minimal_Nbr_Cells_in_Cluster_differential_one_cluster_against_others_Seurat_integration,
Minimal_Nbr_Cells_in_other_Clusters_differential_one_cluster_against_others_Seurat_integration,
logfc_threshold_differential_one_cluster_against_others_Seurat_integration,
Minimal_Percentage_Cells_in_Cluster_differential_one_cluster_against_others_Seurat_integration,
Minimal_Nbr_Cells_in_Cluster_differential_cluster_pairs,
Minimal_Percentage_Cells_in_Cluster_differential_cluster_pairs,
logfc_threshold_differential_cluster_pairs){#accolade ouvrante de la fonction
    
    #chargement des packages
    require(readr)
    require(Seurat)
    require(tidyr)
    require(dplyr)
    require(purrr)
    require(ggplot2)
    require(chameleon)
    
    #conserver le comportement précédent de nest et unnest
    #nest <- nest_legacy
    #unnest <- unnest_legacy
    
    
    #on change la manière d'appeler les dataframes dans les orig.ident des dataframes d'origine projeté dans l'intégration.
    normalized_data_transformed <- lapply(normalized_data,function(x){
        
        newdata <- x
        
        newdata@meta.data$orig.ident <- unlist(lapply(newdata@meta.data$orig.ident,function(z){return(equivalence_nom_sample_nom_figure[which(names(equivalence_nom_sample_nom_figure)==z)])}))
        
        newdata@meta.data$orig.ident <- as.factor(newdata@meta.data$orig.ident)
        
        newdata@meta.data$orig.ident <- droplevels(newdata@meta.data$orig.ident)
        
        return(newdata)
        
    })
    
    #On change également les noms des élément de la list pour que celle-ci puisse être utilisé correctement dans les fonctions qui font appel à elle.
    names(normalized_data_transformed) <- unlist(lapply(names(normalized_data_transformed),function(z){return(equivalence_nom_sample_nom_figure[which(names(equivalence_nom_sample_nom_figure)==z)])}))
    
    
    nom_reduction_transformed <- nom_reduction
    
    names(nom_reduction_transformed) <- unlist(lapply(names(nom_reduction_transformed),function(z){return(equivalence_nom_sample_nom_figure[which(names(equivalence_nom_sample_nom_figure)==z)])}))
    
    
    
    ################code génération répertoire et figures pour mieux comprendre et évaluer les intégrations de Seurat.
    #répertoire résultat
    rootdir <- rootdir_output_function
    
    #Génération du répertoire qui va contenir ls résultats
    dir.create(rootdir)
    
    set.seed(random_seed)
    purrr::pmap(
    .l=list(
    nbdim = as.list(rep(vec_dimension,each=length(vec_resolution))),
    reso = as.list(rep(vec_resolution,times=length(vec_dimension)))
    ),.f = function(nbdim,reso){
        
        #création du répertoire qui va contenir les informations pour une résolution donnée.
        dir.create(paste0(rootdir,"/",as.character(nbdim),"_results"))
        
        ##importation du dataset intégré
        alltogether.comb <- readr::read_rds(integrated_dataset_rds)
        
        
        alltogether.combined <- ScaleData(alltogether.comb, verbose = FALSE)
        alltogether.combined <- Seurat::RunPCA(alltogether.combined, npcs = nbdim, verbose = FALSE)
        
        # t-SNE and Clustering
        alltogether.combined <- Seurat::RunUMAP(alltogether.combined, reduction = "pca", dims = 1:nbdim,n.components = 2L)
        alltogether.combined <- Seurat::FindNeighbors(alltogether.combined, reduction = "pca", dims = 1:nbdim)
        
        #enregistrement du fichier résultat obtenu
        
        readr::write_rds(alltogether.combined,paste0(rootdir,"/",as.character(nbdim),"_results/sobj_integrated_",as.character(nbdim),".rds"))
        
        # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.
        alltogether.combined <- Seurat::FindClusters(alltogether.combined, resolution = reso)
        
        #génération d'un répertoire qui va contenir les figures pour un couple (nbdimensions,resolution) donnée
        
        dir.create(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso))
        
        ############################ A SUPPRIMER ENSUITE ####################################
        
        write_rds(alltogether.combined,paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/sobj_integrated_",as.character(nbdim),"_reso_",reso,".rds"))
        
        ###################################################################################
        
        ### changement des noms des niveaux du facteur orig.ident
        #equivalence_nom_sample_nom_figure
        
        alltogether.combined@meta.data$orig.ident <- unlist(lapply(alltogether.combined@meta.data$orig.ident,function(z){return(equivalence_nom_sample_nom_figure[which(names(equivalence_nom_sample_nom_figure)==z)])}))
        
        alltogether.combined@meta.data$orig.ident <- as.factor(alltogether.combined@meta.data$orig.ident)
        
        alltogether.combined@meta.data$orig.ident <- droplevels(alltogether.combined@meta.data$orig.ident)
        
        #### représentation en umap toutes les données sans discrimination (a compléter demain)
        
        print("génération des représentations umap")
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/all_datasets_integration_clustering_umap_1_2.pdf"),width = 15, height = 15)
        
        print(
        Seurat::DimPlot(alltogether.combined,reduction = "umap",
        pt.size = 1.1,
        cols = chameleon::distinct_colors(length(unique(alltogether.combined@meta.data$seurat_clusters)))[["name"]],
        label = TRUE,
        label.size = 9,
        label.color="red",
        label.box=TRUE,
        dims = c(1,2)
        ))
        
        dev.off()
        
        
        #pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/all_datasets_integration_clustering_umap_1_3.pdf"),width = 15,height = 15)
        
        #print(
        #Seurat::DimPlot(alltogether.combined,reduction = "umap",
        #pt.size = 1.1,
        #cols = chameleon::distinct_colors(length(unique(alltogether.combined@meta.data$seurat_clusters)))[["name"]],
        #label = TRUE,
        #label.size = 9,
        #label.color="red",
        #label.box=TRUE,
        #dims = c(1,3)
        #))
        
        #dev.off()
        
        #pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/all_datasets_integration_clustering_umap_2_3.pdf"),width = 15,height = 15)
        
        #print(
        #Seurat::DimPlot(alltogether.combined,
        #reduction = "umap",
        #pt.size = 1.1,
        #cols = chameleon::distinct_colors(length(unique(alltogether.combined@meta.data$seurat_clusters)))[["name"]],
        #label = TRUE,
        #label.size = 9,
        #label.color="red",
        #label.box=TRUE,
        #dims = c(2,3)
        #))
        
        #dev.off()
        
        
        print("représentation de la compostion des clusters en échantillons et inversement")
        
        #### représentation en umap de l'intégration des données
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/cluster_by_dataset.pdf"),
        width = 35,
        height=30)
        
        print(
        Seurat::DimPlot(alltogether.combined, reduction = "umap", split.by = "orig.ident",label=TRUE,label.size = 7,label.color = "red",cols = chameleon::distinct_colors(length(unique(alltogether.combined@meta.data$seurat_clusters)))[["name"]]
        ))
        
        dev.off()
        
        
        purrr::pmap(.l = list(genelist = gene_pop,nom_pop = as.list(names(gene_pop))),
        .f = purrr::possibly(function(genelist,nom_pop){
            pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/",nom_pop,".pdf"),width = 10,height=10)
            print(
            #il faudrait probablement enlever le min.cutoff
            Seurat::FeaturePlot(alltogether.combined, features = genelist, min.cutoff = "q9",combine = TRUE)
            )
            dev.off()
            
            
            
        },"error : no genes present in the matrix"))
        ## représentation en umap de l'expression d'un certains nombres de gènes
        
        
        print("représentation de différents paramètres importants en Single Cell")
        
        ##  représentation en umap de certaines caractéristiques importantes des données
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/important_parameters.pdf"))
        print(
        Seurat::FeaturePlot(alltogether.combined, features = c("percent_mt"),combine = TRUE)+
        Seurat::DarkTheme()
        )
        
        print(
        Seurat::FeaturePlot(alltogether.combined, features = c("percent_rb"))+
        Seurat::DarkTheme()
        )
        
        print(
        Seurat::FeaturePlot(alltogether.combined, features = c("nFeature_RNA"),combine = TRUE)+
        Seurat::DarkTheme()
        )
        
        print(
        Seurat::FeaturePlot(alltogether.combined, features = c("nCount_RNA"))+
        Seurat::DarkTheme()
        )
        
        
        dev.off()
        
        
        ## exportation d'une table récapitulative de la composition des différents datasets qui ont été intégrés dans les clusters créés
        
    
        print("exportation d'une table récapitulative de la composition des différents datasets qui ont été intégrés dans les clusters créés")
    write.table(df_percentage_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"orig.ident","seurat_clusters"),
        paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/composition_clusters.txt"),quote = FALSE,sep="\t",col.names = TRUE,row.names = FALSE)
        
        
        ## Représentation de la composition des différents datasets qui ont été intégrés dans les clusters créés
        
        
       print("Représentation de la composition des différents datasets qui ont été intégrés dans les clusters créés")
        
        #Définition d un vecteur de couleurs pour qu il y est correspondance entre le code couleurs des heatmaps et de certains barplots
        col_to_use_clusters <- chameleon::distinct_colors(length(unique(alltogether.combined@meta.data$seurat_clusters)))[["name"]]
        
        names(col_to_use_clusters) <- as.character(0:(length(unique(as.character(alltogether.combined@meta.data$seurat_clusters)))-1))
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/composition_clusters.pdf"),
        width = 10,height=10)
        
        #Création d'une graine aléatoire
        set.seed(random_seed)
        
        df_rep_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"orig.ident","seurat_clusters",personnal_colors= TRUE,colchoose = col_to_use_clusters,optionrep = "percentage")
        
        dev.off()
        
        ## Représentation de la composition des différents clusters en goutelettes des différents datasets intégrés (en valeurs brutes : nombre de goutelettes)
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/composition_clusters_valeurs_brutes.pdf"),width = 10,height=10)
        
        df_rep_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"seurat_clusters","orig.ident",optionrep  = "brut")
        
        dev.off()
        
        ## Représentation de la composition des différents clusters en goutelettes des différents datasets intégrés (en pourcentage)
        
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/composition_clusters_percentage.pdf"),width = 10,height=10)
        
        
        df_rep_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"seurat_clusters","orig.ident",optionrep = "percentage")
        
        dev.off()
        
        
        
        #### projection des points de chaque dataset dans le dataset intégré (pts en rouge, les autres points étant en gris)
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/cluster_red_grey.pdf"),width = 30,height=10)
        for(i in 1:(alltogether.combined@meta.data$orig.ident %>% as.factor() %>% nlevels())){
            
            
            
            
            colvec <- character((alltogether.combined@meta.data$orig.ident %>% as.factor() %>% nlevels()))
            for(z in 1:length(colvec)){
                
                if(z==i){
                    
                    colvec[z] <- "red"
                    
                    
                }else{
                    
                    colvec[z] <- "grey"
                    
                }
                
                
                
            }
            
            print(
            Seurat::DimPlot(alltogether.combined,cols = colvec,reduction = "umap",group.by = "orig.ident")
            )
            
            rm(colvec)
            
        }
        dev.off()
        
        #mettre le code ci-dessus également en rouge et blanc.
        
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/cluster_red_white.pdf"),width = 10,height=10)
        for(i in 1:(alltogether.combined@meta.data$orig.ident %>% as.factor() %>% nlevels())){
            
            
            
            
            colvec <- character((alltogether.combined@meta.data$orig.ident %>% as.factor() %>% nlevels()))
            for(z in 1:length(colvec)){
                
                if(z==i){
                    
                    colvec[z] <- "red"
                    
                    
                }else{
                    
                    colvec[z] <- "white"
                    
                }
                
                
                
            }
            
            print(
            DimPlot(alltogether.combined,cols = colvec,reduction = "umap",group.by = "orig.ident")
            )
            
            rm(colvec)
            
        }
        dev.off()
        
        
        
        
        
        #exportation tableau cycle cellulaire
        ##cyclone
        
        write.table(df_percentage_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"seurat_clusters","Cyclone.Phase"),paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/Cyclone_phase_cycle_cellulaire.txt"),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
        
        
        
        ##seurat
        write.table(df_percentage_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"seurat_clusters","Seurat.Phase"),paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/Seurat_phase_cycle_cellulaire.txt"),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
        
        
        #### représentation du cycle cellulaire
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/phase_cycle_cellulaire.pdf"),width = 10,height=10)
        ##représentation en umap (phase Cyclone)
        print(
        
        Seurat::DimPlot(alltogether.combined,reduction="umap",group.by = "Cyclone.Phase")
        
        )
        ##représentation en umap (phase Seurat)
        print(
        
        Seurat::DimPlot(alltogether.combined,reduction="umap",group.by = "Seurat.Phase")
        
        )
        
        set.seed(random_seed)
        ## représentation cycle cellulaire
        ##représentation histogramme (stade dans le cycle cellulaire prédit par Cyclone des cellules de chaque cluster)
        
        
        
        df_rep_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"seurat_clusters","Cyclone.Phase",personnal_colors= TRUE, colchoose = c("G1"="blue","G2M"="green","S" = "red"),optionrep = "percentage")
        
        ##représentation histogramme (stade dans le cycle cellulaire prédit par Seurat des cellules de chaque cluster)
        
        set.seed(random_seed)
        
        
        df_rep_of_each_level_a_factor_in_each_level_another_factor(alltogether.combined@meta.data,"seurat_clusters","Seurat.Phase",personnal_colors= TRUE, colchoose = c("G1"="blue","G2M"="green","S" = "red"),optionrep = "percentage")
        
        dev.off()
        
        
        
        #On visualise où les goutelettes se trouvent dans les datasets d'origine
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/droplets_in_original_datasets.pdf"),width = 10,height=10)
        
        plotclusterpoint_in_original_dataset(alltogether.combined,normalized_data_transformed,nom_reduction_transformed)
        
        
        dev.off()
        
        #annotation du jeu de données avec SingleR
        print("début Analyses SingleR")
        
        representation_annotation_SingleR_Seurat_Integration(normalized_data_transformed,alltogether.combined,ref_SingleR,vector_col_annotation = vector_col_annota,paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/annotation_SingleR_umap.pdf"),height_pdf = 15, width_pdf = 15,pt_size_dimplot = pt_size_dimplot_SingleR,path_rds_file = paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/annotation_SingleR_metadata.rds"))
        
        #lecture du fichier qui contient les stats SingleR
        stat_SingleR <- read_rds(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/annotation_SingleR_metadata.rds"))
        
        pdf(paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/SingleR_composition_type_cellulaire_par_cluster.pdf"),height = 15 , width=15)
        
        df_rep_of_each_level_a_factor_in_each_level_another_factor(stat_SingleR,"seurat_clusters","all_single_R_annotation",personnal_colors= TRUE, colchoose = vector_col_annota,"percentage")
        
        dev.off()
        
        rm(stat_SingleR)
        
        print("fin analyses SingleR")
        
       # print("analyses différentielles Un cluster vs tous les autres")
        
      #  differential_one_cluster_against_others_Seurat_integration(normalized_individual_datasets = normalized_data_transformed,Seurat_Integration = alltogether.combined,
       # Minimal_Nbr_Cells_in_Cluster = Minimal_Nbr_Cells_in_Cluster_differential_one_cluster_against_others_Seurat_integration,
       # Minimal_Nbr_Cells_in_other_Clusters = Minimal_Nbr_Cells_in_other_Clusters_differential_one_cluster_against_others_Seurat_integration,
       # logfc_threshold = logfc_threshold_differential_one_cluster_against_others_Seurat_integration,
       # path_results = paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/differential_analysis_one_cluster_against_all_others"),
       # Minimal_Percentage_Cells_in_Cluster = Minimal_Percentage_Cells_in_Cluster_differential_one_cluster_against_others_Seurat_integration)
        
      #  print("analyses différentielles : toutes les comparaisons de 2 clusters possibles")
        
        
       # differential_cluster_pairs(normalized_individual_datasets = normalized_data_transformed,Seurat_Integration = alltogether.combined,Minimal_Nbr_Cells_in_Cluster = #Minimal_Nbr_Cells_in_Cluster_differential_cluster_pairs ,Minimal_Percentage_Cells_in_Cluster = Minimal_Percentage_Cells_in_Cluster_differential_cluster_pairs,logfc_threshold = #logfc_threshold_differential_cluster_pairs,path_results = paste0(rootdir,"/",as.character(nbdim),"_results/nbdim_",nbdim,"_reso_",reso,"/differential_analysis_comparisons_two_clusters"))
        
        ##suppression de alltogether.comb (par prudence mais ce n'est peut être pas nécessaire)
        rm(alltogether.comb)
        #suppression de l'objet Seurat avec les datasets intégrés pour les valeurs de résolutions et de nombre de dimensions fixés.
        rm(alltogether.combined)
    })
    
    
}#accolade fermante de la fonction

