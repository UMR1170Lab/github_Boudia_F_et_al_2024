
#vidange de l environnement
rm(list=ls())

require(purrr)
require(DESeq2)
require(GenomicFeatures)
require(tximport)

#Loading of some functions that are use in this analysis
source("/path/to/LibrarySizeBarplot.R")
source("/path/to/PCATSNErep.R")
source("/path/to/DESeq2_resultdup_R4.R")

## paths fichiers pour le traitemen des données / chaînes de caractères spécifiques pour certaines commandes
PathFileSFQuelconque <- "/path/to/the/directory/sf_files/any_sf_files.sf.sf"
PathFileGTF <- "/path/to/the/directory/gtf_file/Homo_sapiens.GRCh38.97/Homo_sapiens.GRCh38.97.gtf"# peut-être à changer
PathInternetFileGTF <- "ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/"
path_sf_file <- "/path/to/the/directory/sf_files/"
StudiedOrganism <- "Homo sapiens"
AlignmentTool <- "salmon"

#importation du fichier GTF dans R
gtf_file <- rtracklayer::import(paste0(PathFileGTF))
# class(gtf_file)
# colnames(gtf_file@elementMetadata)


#génération d un dataframe permettant de passer de l identifiant ensembl au gene symbol facilement.
gtf_file@elementMetadata[grep("CTNNB1",gtf_file@elementMetadata$gene_name),]
genes_conv <- gtf_file@elementMetadata
genes_conv <- as.data.frame(genes_conv)
colnames(genes_conv)




# Conversion du GTF dans un format qui va pouvoir être utilisé par txiimport et probablement DeSeq2 aussi.
tx2gtf <- makeTxDbFromGFF(file= paste0(PathFileGTF),format = "auto",dataSource = paste0(PathInternetFileGTF),organism = paste0(StudiedOrganism))
tx2gtf <- as.list(tx2gtf)
tx2gene <- data.frame(tx_id = tx2gtf$transcripts$tx_name,gene_id = tx2gtf$genes$gene_id)
# head(tx2gene)
dim(tx2gene)
tx2gene$tx_id <- as.character(tx2gene$tx_id)
tx2gene$gene_id <-as.character(tx2gene$gene_id)


#section permettant de spécifier les paths des fichiers de quantification Salmon (ainsi que la condition à laquelle chacun est associé) en vue de leur utilisation dans R.
files <- unname(sapply(dir(path_sf_file)[1:length(dir(path_sf_file))],FUN = function(x) paste0(path_sf_file,x)))# vecteur qui contient le nom des fichiers
names(files) <- dir(path_sf_file)[1:length(dir(path_sf_file))]#on affecte un nom à chaque élément du vecteur files qui va correspondre aux noms des colonnes de tximport
names(files) <- unlist(sapply(names(files),FUN = function(x){strsplit(x,".sf","")}))

names(files)
files

#génération d un dataframe de metadata spécifiant pour chaque échantillon à quelle condition il est associé.
metadata <- data.frame(condition = c(
rep("MO7e_shDLX3",times=3),
rep("MO7e_shRen",times=3),
rep("WSU_shRen",times=3),
rep("Wsu_shDLX3",times=3)
))
#les noms de lignes de ce dataframe vont correspondre aux noms du vecteur files
rownames(metadata) <- names(files)

#On peut ensuite employer la fonction tximport pour importer les fichiers de quantification dans R
txi <- tximport(files,type = "salmon",txOut = FALSE,tx2gene=tx2gene,ignoreTxVersion=TRUE)

#affichage de la matrice de counts de txi
dim(txi$counts)

#spécification des couleurs associées à chaque condition.
colorsample <- c("#3E67C5","#9C3FC6","#9B5188","#C6314D","#30E36E","#4DD6FB","#51DDB6")
LibrarySizeBarplot(Data = txi$counts,Metadata = metadata,ColMetadata = condition,Couleur = colorsample,xlabtitle = "samples",ylabtitle = "taille des librairies",TitreFigure = "Taille des librairies",TransLog = FALSE)
dds <- DESeqDataSetFromTximport(txi, metadata, ~condition)
## sélection des gènes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
## Normalisation de la matrice
### estimation des facteurs de taille
dds <- estimateSizeFactors(dds)
#affichage des facteurs de taille
sizeFactors(dds)
#génération d une première matrice normalisée.
normalized_counts <- counts(dds, normalized=TRUE)

#affichage des dimensions de la matrice normalisée.
dim(normalized_counts)
head(normalized_counts)
#Génération fichiers pour GSEA
# Génération des fichiers pour faire des GSEA
## Génération du dataframe pour faire les GSEA

######################################################## première façon de générer une matrice utilisable avec GSEA ##################################
normalized_counts_for_GSEA <- log2(normalized_counts+1)
normalized_counts_for_GSEA <- conv_gene_tab(normalized_counts_for_GSEA,genesconvTab = genes_conv)
head(normalized_counts_for_GSEA)
normalized_counts_for_GSEA <- data.frame(NOM = normalized_counts_for_GSEA$symbol_gene,DESCRIPTION = rep("na",times = nrow(normalized_counts_for_GSEA)),normalized_counts_for_GSEA[,-which(colnames(normalized_counts_for_GSEA)=="symbol_gene")])
colnames(normalized_counts_for_GSEA) <- c("NOM","DESCRIPTION",colnames(normalized_counts))
head(normalized_counts_for_GSEA)
normalized_counts_for_GSEA <- log2(normalized_counts+1)
normalized_counts_for_GSEA <- conv_gene_tab(normalized_counts_for_GSEA,genesconvTab = genes_conv)
head(normalized_counts_for_GSEA)
normalized_counts_for_GSEA <- data.frame(NOM = normalized_counts_for_GSEA$symbol_gene,DESCRIPTION = rep("na",times = nrow(normalized_counts_for_GSEA)),normalized_counts_for_GSEA[,-which(colnames(normalized_counts_for_GSEA)=="symbol_gene")])
colnames(normalized_counts_for_GSEA) <- c("NOM","DESCRIPTION",colnames(normalized_counts))
write.table(normalized_counts_for_GSEA,"/path/to/GSEA_files/normalized_counts_for_GSEA_log2_transformed.txt",row.names=FALSE,quote = FALSE,sep ="\t",col.names = TRUE)
write.table(chipfile,file = paste0("/path/to/GSEA_files/normalized_counts_for_GSEA_log2_transformed.chip"),quote=FALSE,row.names=FALSE,col.names = TRUE,sep="\t")

chipfile <- data.frame("Probe Set ID"=unique(genes_conv$gene_name),
"Gene Symbol" = unique(genes_conv$gene_name),
"Gene Title"=rep("NA",times = length(unique(genes_conv$gene_name))))
 
colnames(chipfile) <- c("Probe Set ID","Gene Symbol","Gene Title")
 
write.table(chipfile,file = paste0("/path/to/GSEA_files/normalized_counts_for_GSEA_log2_transformed.chip"),quote=FALSE,row.names=FALSE,col.names = TRUE,sep="\t")


############################################### partie nouvelle version normalisation matrices et génération des fichiers pour GSEA ##############

################ différents types de normalisation sont testés #############


#génération d une list contenant les différents types de normalisation possible.
matrix_normalization_different_types <- list(counts_dds_normalized = counts(dds, normalized=TRUE),
                                             vst_dds_normalized = vst(dds,blind = FALSE),
                                             rlog_dds_normalized = rlog(dds, blind=FALSE))



#représentation des données quand elles sont normalisées avec la méthode vst.

pdf("/path/to/output_files_directory/representation_pca_tsne_vst.pdf",width = 10,height = 10)

PCATSNErep(Data = as.matrix(t(assay(matrix_normalization_different_types$vst_dds_normalized))),
Metadata = metadata,ColMetadata = condition,
PCA_Arg = list(scaling = TRUE,ndim = 3, DimX = "Dim.1",DimY = "Dim.2",title = "PCA"),
TSNE_arg = list(ndim = 3,perp = 3,nb_itermax =4000),
ColorPoint = colorsample)

dev.off()

#représentation des données quand elles sont normalisées avec la méthode rlog.

pdf("/path/to/output_files_directory/representation_pca_tsne_rlog.pdf")

PCATSNErep(Data = as.matrix(t(assay(matrix_normalization_different_types$rlog_dds_normalized))),
Metadata = metadata,ColMetadata = condition,
PCA_Arg = list(scaling = TRUE,ndim = 3, DimX = "Dim.1",DimY = "Dim.2",title = "PCA"),
TSNE_arg = list(ndim = 3,perp = 3,nb_itermax =4000),
ColorPoint = colorsample)


dev.off()


########## on s occupe de générer une matrice normalisée avec la méthode vst 

vst <- vst(dds, blind=FALSE) 

normalized_counts_vst_for_GSEA <- assay(vst)

PathExpressionFileGSEA_vst <- "/path/to/GSEA_files_directory/normalized_expression_vst_for_GSEA.txt"
PathMetadataGSEA_vst <- "/path/to/GSEA_files_directory/metadata_for_gsea.cls"
PathChipFileGSEA_vst <- "/path/to/GSEA_files_directory/chipfile_for_gsea.chip"



normalized_counts_vst_for_GSEA <- conv_gene_tab(normalized_counts_vst_for_GSEA,genesconvTab = genes_conv)
normalized_counts_vst_for_GSEA <- data.frame(NOM = normalized_counts_vst_for_GSEA$symbol_gene,DESCRIPTION = rep("na",times = nrow(normalized_counts_vst_for_GSEA)),normalized_counts_vst_for_GSEA[,-which(colnames(normalized_counts_vst_for_GSEA)=="symbol_gene")])
colnames(normalized_counts_vst_for_GSEA) <- c("NOM","DESCRIPTION",colnames(assay(vst)))
 write.table(normalized_counts_vst_for_GSEA,paste0(PathExpressionFileGSEA_vst),row.names=FALSE,quote = FALSE,sep ="\t",col.names = TRUE)
capture.output(
    
 # ############## Chaîne de caractères pour faire le fichier cls ##############################
  cat(cat(c(length(metadata$condition),nlevels(as.factor(metadata$condition)),"1"),"\n#",c(unique(as.character(metadata$condition))),"\n"),cat(as.character(metadata$condition)),sep=""),file = paste0(PathMetadataGSEA_vst)
 
 
 
  )
chipfile <- data.frame("Probe Set ID"=unique(genes_conv$gene_name),
"Gene Symbol" = unique(genes_conv$gene_name),
"Gene Title"=rep("NA",times = length(unique(genes_conv$gene_name))))

colnames(chipfile) <- c("Probe Set ID","Gene Symbol","Gene Title")

write.table(chipfile,file = paste0(PathChipFileGSEA_vst),quote=FALSE,row.names=FALSE,col.names = TRUE,sep="\t")





################################################## partie analyses différentielles #####################################################

#on lance les analyses différentielles.
difftest <- DEseq2_resultdup_R4(txi,metadataTab = metadata,threshold_keep_gene=10)

difftest_avecShrinkage <- lapply(difftest,function(z){return(z[["avecShrinkage"]])})

head(difftest_avecShrinkage$Wsu_shDLX3_vs_WSU_shRen)

head(difftest_avecShrinkage$Wsu_shDLX3_vs_WSU_shRen,n=100)

#exportation des analyses différentielles (mode avec shrinkage et mode sans shrinkage). Ici au format rds.
write_rds(difftest,"/WORKDIR/media/umr1170/saturne/RNAseq_fabien_sh_DLX3_human_cell_lineages/file_differentiel_rds/file_all_differentiel.rds")
#exportation du mode avec Shrinkage.
write_rds(difftest_avecShrinkage,"/WORKDIR/media/umr1170/saturne/RNAseq_fabien_sh_DLX3_human_cell_lineages/file_differentiel_rds/file_only_withShrinkage_differentiel.rds")

#exportation au format xls des résultats d analyses différentielles.
WriteXLS::WriteXLS(difftest_avecShrinkage,SheetNames = names(difftest_avecShrinkage),row.names=TRUE,col.names=TRUE,ExcelFileName="/WORKDIR/media/umr1170/saturne/RNAseq_fabien_sh_DLX3_human_cell_lineages/file_differentiel_excel/differentiel_avec_shrinkage.xls")


