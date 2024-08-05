#!/bin/bash


########################################### Objectif du script ###############################################################
#l'objectif de ce script est de faire une annotation des fichiers bed issues de macs2 callpeak

###Variables dont l'utilisateur doit définir les valeurs 


#Path contenant les fichiers txt qui doivent être annotés (proviennent des fichiers bed issus de diffbind)
pathtxtfile="/path/to/analysis_all_samples/deseq2/vs/FDR_001_sorted/CB_CTRLvsCB_DLX3"


#path où les fichiers d'annotations vont être générés.
pathoutput="/path/to/analysis_all_samples/deseq2/vs/annotations_motifs_peaks_up_down_peaks_diff_LFC_sup1_FDR_001"

#nombre de threads
nbThreads=15

############################

#déplacements dans le répertoire où se trouvent les fichiers qui vont être annotés
cd $pathtxtfile

for i in `ls *.bed`
do

#extraction du prefixe du fichier txt
prefixfile=$(echo $i | sed 's/.bed//')

echo "le prefix de ce fichier txt est ${prefixfile}"

#génération d'un dossier qui contiendra les résultats des analyses associées à ce fichier txt

mkdir ${pathoutput}/${prefixfile}_dossier

#génération d'un dossier qui contiendra les résultats des annotations de peaks associées à ce fichier txt.
mkdir ${pathoutput}/${prefixfile}_dossier/dossier_annotations

#génération d'un dossier qui contiendra les résultats des analyses d'enrichissement GenomeOntology

mkdir ${pathoutput}/${prefixfile}_dossier/dossier_annotations/GeneOntology_enrichresult

#génération d'un dossier qui contiendra les résultats des analyses d'enrichissement GO

mkdir ${pathoutput}/${prefixfile}_dossier/dossier_annotations/go_enrichresult

#génération d'un dossier qui contiendra les résultats des analyses de motifs associées à ce fichier txt.
mkdir ${pathoutput}/${prefixfile}_dossier/dossier_motif


echo "Analyse de motifs"

#analyses de motifs
findMotifsGenome.pl $i hg38 ${pathoutput}/${prefixfile}_dossier/dossier_motif -p $nbThreads 
#-size 75
#annotations des peaks

echo "Annotations"
#-genomeOntology ${pathoutput}/${prefixfile}_dossier/dossier_annotations/GeneOntology_enrichresult
#-go ${pathoutput}/${prefixfile}_dossier/dossier_annotations/go_enrichresult
#-m ${pathoutput}/${prefixfile}_dossier/dossier_motif/*.all.motifs
#-annStats ${pathoutput}/${prefixfile}_dossier/dossier_annotations/${prefixfile}_Stats.tsv
#/home/umr1170/anaconda2/envs/chipseq_v1/share/homer/data/genomes/hg38
annotatePeaks.pl ${i} hg38  > ${pathoutput}/${prefixfile}_dossier/dossier_annotations/${prefixfile}.tsv

done

#-genomeOntology ${pathoutput}/${prefixfile}_dossier/dossier_annotations/GeneOntology_enrichresult \
#-go ${pathoutput}/${prefixfile}_dossier/dossier_annotations/go_enrichresult -annStats ${pathoutput}/${prefixfile}_dossier/dossier_annotations/${prefixfile}_Stats.tsv



