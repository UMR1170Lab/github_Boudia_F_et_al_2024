require(ggplot2)
require(ggpie)
require(readr)
show_samples_pca_coordinate_file <- function(path_tsv_file_pca_coordinates){#accolade ouvrante de la fonction
# importation du fichier tsv contenant les coordonnées pca associés à chaque échantillon.
df <- read.table(path_tsv_file_pca_coordinates,skip= 11,h = T)
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _
df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0
(z[-length(z)],collapse="_")}))
print("affichage des premières lignes du dataframe")
print(head(df))
print("affichage des conditions associés à des points (des échantillons) dans la pca")
print(unique(df$condition))
#accolade fermante de la fonction
}
customized_representation_pca_from_nfcore_atacseq_pipeline <- function(path_tsv_file_pca_coordinates,colors_condition_in_pca,path_output_pdf,width_pdf =10,height_pdf = 10,size_label_ggrepel = 2,max_iteration_ggrepel = 6000){#accolade ouvrante de la fonction
#chargement des packages nécessaires pour le fonctionnement de la fonction
require(ggplot2)
require(ggrepel)
#importation du fichier tsv contenant les coordonnées PCA des échantillons analysés
df <- read.table(path_tsv_file_pca_coordinates,skip=11,h=T)
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _
df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0(z[-length(z)],collapse="_")}))
#pourcentage de variance associé à l axe 1
percent_variance_PC1 <- strsplit(colnames(df)[2],"[..]")[[1]][3]
#pourcentage de variance associé à l axe 2
percent_variance_PC2 <- strsplit(colnames(df)[3],"[..]")[[1]][3]
#on change ensuite les noms des colonnes 2 et 3 (celles dans lesquelles figurent les coordonnées des échantillons pour les axes 1 et 2)
colnames(df)[2] <- "PC1"
colnames(df)[3] <- "PC2"
#affichage des premières lignes du dataframe
print(head(df))
#exportation au format pdf de la figure ggplot (représentation pca avec annotation des échantillons)
pdf(path_output_pdf,width = width_pdf,height = height_pdf)
print(
ggplot2::ggplot(df,aes(x = PC1,y = PC2, color=condition,label = sample))+geom_po
int()+ scale_color_manual(values= colors_condition_in_pca)+
ggrepel::geom_text_repel(size=2,max.iter = max_iteration_ggrepel,label.size =siz
e_label_ggrepel)+
xlab(paste0("PC1 : Variance = ",percent_variance_PC1, " %"))+
ylab(paste0("PC2 : Variance = ",percent_variance_PC2, " %"))+
theme_classic()
)
dev.off()
#accolade fermante de la fonction
}
show_samples_pca_coordinate_file("/path/to/consensus_peaks.mLb.clN.pca.vals_mqc.tsv")
customized_representation_pca_from_nfcore_atacseq_pipeline(path_tsv_file_pca_coordinates = "/path/to/coordinate/consensus_peaks.mLb.clN.pca.vals_mqc.tsv",
colors_condition_in_pca = c("CB1_CTRL" = "black","CB2_CTRL" = "black","CB1_DLX3"= "green4","CB2_DLX3" = "palegreen2"),
path_output_pdf = "/path/to/output/directory/representation_PCA_dim12_v1.pdf",
width_pdf =10,
height_pdf = 10,
size_label_ggrepel = 2,
max_iteration_ggrepel = 6000)

customized_representation_pca_from_nfcore_atacseq_pipeline <- function(path_tsv_file_pca_coordinates,colors_condition_in_pca,path_output_pdf,
width_pdf =10,height_pdf = 10,
size_label_ggrepel = 2,
max_iteration_ggrepel = 6000){#accolade ouvrante de la fonction
#chargement des packages nécessaires pour le fonctionnement de la fonction
require(ggplot2)
require(ggrepel)
#importation du fichier tsv contenant les coordonnées PCA des échantillons analysés
df <- read.table(path_tsv_file_pca_coordinates,skip=11,h=T,sep = "\t")
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _
df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0(z[-length(z)],collapse="_")}))
#pourcentage de variance associé à l axe 1
percent_variance_PC1 <- strsplit(colnames(df)[2],"[..]")[[1]][3]
#pourcentage de variance associé à l axe 2
percent_variance_PC2 <- strsplit(colnames(df)[3],"[..]")[[1]][3]
#on change ensuite les noms des colonnes 2 et 3 (celles dans lesquelles figurent les coordonnées des échantillons pour les axes 1 et 2)
colnames(df)[2] <- "PC1"
colnames(df)[3] <- "PC2"
#affichage des premières lignes du dataframe
print(head(df))
#exportation au format pdf de la figure ggplot (représentation pca avec annotation des échantillons)
pdf(path_output_pdf,width = width_pdf,height = height_pdf)
print(
ggplot2::ggplot(df,aes(x = PC1,y = PC2, color=condition,label = sample))+
geom_point()+
scale_color_manual(values= colors_condition_in_pca)+
ggrepel::geom_text_repel(size=2,max.iter = max_iteration_ggrepel,label.size =size_label_ggrepel)+
xlab(paste0("PC1 : Variance = ",percent_variance_PC1, " %"))+
ylab(paste0("PC2 : Variance = ",percent_variance_PC2, " %"))+
theme_classic()
)
dev.off()
#accolade fermante de la fonction
}
customized_representation_pca_from_nfcore_atacseq_pipeline(path_tsv_file_pca_coordinates = "/path/to/pca/coordinates/consensus_peaks.mLb.clN.pca.vals_mqc.tsv",

colors_condition_in_pca = c("CB1_CTRL" = "black","CB2_CTRL" = "black","CB1_DLX3" = "green4","CB2_DLX3" = "palegreen2"),
path_output_pdf = "/path/to/output/directory/representation_PCA_dim12_v1.pdf",
width_pdf =10,
height_pdf = 10,
size_label_ggrepel = 2,
max_iteration_ggrepel = 6000)

customized_representation_pca_from_nfcore_atacseq_pipeline <- function(path_tsv_file_pca_coordinates,colors_condition_in_pca,path_output_pdf,width_pdf =10,height_pdf = 10,size_label_ggrepel = 2,max_iteration_ggrepel = 6000){#accolade ouvrante de la fonction
#chargement des packages nécessaires pour le fonctionnement de la fonction
require(ggplot2)
require(ggrepel)
#importation du fichier tsv contenant les coordonnées PCA des échantillons analysés
df <- read.table(path_tsv_file_pca_coordinates,skip=11,h=T,sep = "\t")
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _
df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0(z[-length(z)],collapse="_")}))
#pourcentage de variance associé à l axe 1
percent_variance_PC1 <- strsplit(colnames(df)[2],"[..]")[[1]][3]
#pourcentage de variance associé à l axe 2
percent_variance_PC2 <- strsplit(colnames(df)[3],"[..]")[[1]][3]
#on change ensuite les noms des colonnes 2 et 3 (celles dans lesquelles figurent les coordonnées des échantillons pour les axes 1 et 2)
colnames(df)[2] <- "PC1"
colnames(df)[3] <- "PC2"
#affichage des premières lignes du dataframe
print(head(df))
#exportation au format pdf de la figure ggplot (représentation pca avec annotation des échantillons)
pdf(path_output_pdf,width = width_pdf,height = height_pdf)
print(
ggplot2::ggplot(df,aes(x = PC1,y = PC2, color=condition,label = sample))+geom_point(size = 2)+ scale_color_manual(values= colors_condition_in_pca)+
ggrepel::geom_text_repel(size=2,max.iter = max_iteration_ggrepel,label.size =size_label_ggrepel)+
xlab(paste0("PC1 : Variance = ",percent_variance_PC1, " %"))+
ylab(paste0("PC2 : Variance = ",percent_variance_PC2, " %"))+
theme_classic()
)
dev.off()
#accolade fermante de la fonction
}
customized_representation_pca_from_nfcore_atacseq_pipeline(path_tsv_file_pca_coordinates = "/path/to/consensus_peaks.mLb.clN.pca.vals_mqc.tsv",
colors_condition_in_pca = c("CB1_CTRL" = "black","CB2_CTRL" = "black","CB1_DLX3" = "green4","CB2_DLX3" = "palegreen2"),
path_output_pdf = "/path/to/output/directory/representation_PCA_dim12_v1.pdf",
width_pdf =10,
height_pdf = 10,
size_label_ggrepel = 2,
max_iteration_ggrepel = 6000)

customized_representation_pca_from_nfcore_atacseq_pipeline <- function(path_tsv_file_pca_coordinates,colors_condition_in_pca,path_output_pdf,width_pdf =10,height_pdf = 10,size_label_ggrepel = 2,max_iteration_ggrepel = 6000){#accolade ouvrante de la fonction

#chargement des packages nécessaires pour le fonctionnement de la fonction
require(ggplot2)
require(ggrepel)
#importation du fichier tsv contenant les coordonnées PCA des échantillons analysés
df <- read.table(path_tsv_file_pca_coordinates,skip=11,h=T,sep = "\t")
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _
df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0(z[-length(z)],collapse="_")}))
#pourcentage de variance associé à l axe 1
percent_variance_PC1 <- strsplit(colnames(df)[2],"[..]")[[1]][3]
#pourcentage de variance associé à l axe 2
percent_variance_PC2 <- strsplit(colnames(df)[3],"[..]")[[1]][3]
#on change ensuite les noms des colonnes 2 et 3 (celles dans lesquelles figurent les coordonnées des échantillons pour les axes 1 et 2)
colnames(df)[2] <- "PC1"
colnames(df)[3] <- "PC2"
#affichage des premières lignes du dataframe
print(head(df))
#exportation au format pdf de la figure ggplot (représentation pca avec annotation des échantillons)

pdf(path_output_pdf,width = width_pdf,height = height_pdf)
print(
ggplot2::ggplot(df,aes(x = PC1,y = PC2, color=condition,label = sample))+geom_point(size = 4)+ scale_color_manual(values= colors_condition_in_pca)+
ggrepel::geom_text_repel(size=2,max.iter = max_iteration_ggrepel,label.size =siz
e_label_ggrepel)+
xlab(paste0("PC1 : Variance = ",percent_variance_PC1, " %"))+
ylab(paste0("PC2 : Variance = ",percent_variance_PC2, " %"))+
theme_classic()
)
dev.off()
#accolade fermante de la fonction
}
customized_representation_pca_from_nfcore_atacseq_pipeline(path_tsv_file_pca_coordinates = "/path/to/consensus_peaks.mLb.clN.pca.vals_mqc.tsv",
colors_condition_in_pca = c("CB1_CTRL" = "black","CB2_CTRL" = "black","CB1_DLX3" = "green4","CB2_DLX3" = "palegreen2"),
path_output_pdf = "/path/to/output/directory/representation_PCA_dim12_v1.pdf",
width_pdf =5,
height_pdf = 5,
size_label_ggrepel = 2,
max_iteration_ggrepel = 6000)
customized_representation_pca_from_nfcore_atacseq_pipeline <- function(path_tsv_file_pca_coordinates,colors_condition_in_pca,path_output_pdf,width_pdf =10,height_pdf = 10,size_label_ggrepel = 2,max_iteration_ggrepel = 6000){#accolade ouvrante de la fonction
#chargement des packages nécessaires pour le fonctionnement de la fonction
require(ggplot2)
require(ggrepel)
#importation du fichier tsv contenant les coordonnées PCA des échantillons analysés
df <- read.table(path_tsv_file_pca_coordinates,skip=11,h=T,sep = "\t")
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _

df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0(z[-length(z)],collapse="_")}))

#pourcentage de variance associé à l axe 1
percent_variance_PC1 <- strsplit(colnames(df)[2],"[..]")[[1]][3]
#pourcentage de variance associé à l axe 2
percent_variance_PC2 <- strsplit(colnames(df)[3],"[..]")[[1]][3]
#on change ensuite les noms des colonnes 2 et 3 (celles dans lesquelles figurent les coordonnées des échantillons pour les axes 1 et 2)
colnames(df)[2] <- "PC1"
colnames(df)[3] <- "PC2"
#affichage des premières lignes du dataframe
print(head(df))
#exportation au format pdf de la figure ggplot (représentation pca avec annotation des échantillons)
pdf(path_output_pdf,width = width_pdf,height = height_pdf)
print(
ggplot2::ggplot(df,aes(x = PC1,y = PC2, color=condition,label = sample))+geom_point(size = 2.5)+
scale_color_manual(values= colors_condition_in_pca)+
ggrepel::geom_text_repel(size=2,max.iter = max_iteration_ggrepel,label.size =size_label_ggrepel)+
xlab(paste0("PC1 : Variance = ",percent_variance_PC1, " %"))+
ylab(paste0("PC2 : Variance = ",percent_variance_PC2, " %"))+
theme_classic()
)
dev.off()
#accolade fermante de la fonction
}

customized_representation_pca_from_nfcore_atacseq_pipeline(path_tsv_file_pca_coordinates = "/path/to/consensus_peaks.mLb.clN.pca.vals_mqc.tsv",
colors_condition_in_pca = c("CB1_CTRL" = "black","CB2_CTRL" = "black","CB1_DLX3" = "green4","CB2_DLX3" = "palegreen2"),
path_output_pdf = "/path/to/output/directory/representation_PCA_dim12_v1.pdf",
width_pdf =5,
height_pdf = 5,
size_label_ggrepel = 2,
max_iteration_ggrepel = 6000)
customized_representation_pca_from_nfcore_atacseq_pipeline <- function(path_tsv_file_pca_coordinates,colors_condition_in_pca,path_output_pdf,width_pdf =10,height_pdf = 10,max_iteration_ggrepel = 6000){#accolade ouvrante de la fonction
#chargement des packages nécessaires pour le fonctionnement de la fonction
require(ggplot2)
require(ggrepel)
#importation du fichier tsv contenant les coordonnées PCA des échantillons analysés
df <- read.table(path_tsv_file_pca_coordinates,skip=11,h=T,sep = "\t")
#génération d une colonne condition qui donne pour chaque échantillon, la condition à laquelle il est associé. Il faut que les espaces dans les noms des conditions soient obligatoirement associés à des _
df$condition <- unlist(lapply(strsplit(unique(df$sample),"_"),function(z){paste0(z[-length(z)],collapse="_")}))
#pourcentage de variance associé à l axe 1
percent_variance_PC1 <- strsplit(colnames(df)[2],"[..]")[[1]][3]
#pourcentage de variance associé à l axe 2
percent_variance_PC2 <- strsplit(colnames(df)[3],"[..]")[[1]][3]
#on change ensuite les noms des colonnes 2 et 3 (celles dans lesquelles figurent les coordonnées des échantillons pour les axes 1 et 2)
colnames(df)[2] <- "PC1"
colnames(df)[3] <- "PC2"
#affichage des premières lignes du dataframe
print(head(df))
#exportation au format pdf de la figure ggplot (représentation pca avec annotation des échantillons)
pdf(path_output_pdf,width = width_pdf,height = height_pdf)
print(
ggplot2::ggplot(df,aes(x = PC1,y = PC2, color=condition,label = sample))+geom_point(size = 2.5)+ scale_color_manual(values= colors_condition_in_pca)+
xlab(paste0("PC1 : Variance = ",percent_variance_PC1, " %"))+
ylab(paste0("PC2 : Variance = ",percent_variance_PC2, " %"))+
theme_classic()
)
dev.off()
#accolade fermante de la fonction
}

customized_representation_pca_from_nfcore_atacseq_pipeline(path_tsv_file_pca_coordinates = "/path/to/consensus_peaks.mLb.clN.pca.vals_mqc.tsv",
colors_condition_in_pca = c("CB1_CTRL" = "black","CB2_CTRL" = "grey","CB1_DLX3" = "green4","CB2_DLX3" = "palegreen2"),
path_output_pdf = "/path/to/output/directory/representation_PCA_dim12_v1.pdf",
width_pdf =5,
height_pdf = 5,
max_iteration_ggrepel = 6000)

history(max.show=Inf)


