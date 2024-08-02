#!/bin/bash
#SBATCH --job-name=salmon               # Nom du Job
#SBATCH --ntasks=1                       # Nombre de Tasks : 1
#SBATCH --cpus-per-task=14                # Allocation de 12 CPU par Task
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --partition=mediumq
#SBATCH --output=rapport_quantification_salmon.out

#### Objectif 
# effectuer la quantification des fastqs trimmés avec salmon 0.14.2

# Variables à modifier si nécessaire.

## Variable path pour l'activation de conda avec source

path_conda="/path/to/conda.sh"


## nom de l'environnement conda qui contient la version de Salmon que l on va utiliser.

conda_env="/path/to/conda_env/dupsalmonlatest"


## path qui contient les fastqs qui vont être quantifier.

path_fastq="/path/to/input_fastqs"


## extension des fichiers R1

extR1="_val_1.fq.gz"
extR2="_val_2.fq.gz"


##dossier qui va contenir les résultats de la quantifivcation avec Salmon (quantification_fastq_trimmed_with_trimmgalore)
pathresults="/mnt/beegfs/userdata/e_robert/DLX3_sh_RNAseq/output_quantification_salmon"


## Nombre de threads qui vont être utilisés

nbThread=13

############## paramètres pour l'utilisation de Salmon #################################

salmon_index="/path/to/decoy_transcriptome_salmon/human_release_97_ensembl/human_97_ensembl_salmon_0142_decoys_index_k_eg_31"


######################################### code qui va être utilisé #######################################

##activation de conda

source $path_conda

##activation de l'environnement conda contenant salmon.

conda activate $conda_env

## déplacement dans le répertoire contenant les fastqs

cd $path_fastq

## Création d'un array contenant les préfixes des fichiers fastqs qui von être quantifié (un préfixe pour une paire)

arrayfastq=(`find ./ -name "*$extR1" -exec basename {}  \; | sed "s/$extR1//"`)

## variable qui donne la valeur associée à la dernière itération de la boucle qui va être utilisé pour faire la quantification de chacune des paires de fastqs

LengthfileCdOneminusOne=$((${#arrayfastq[@]}-1))


## boucle pour la génération des fichiers de quantification

for i in `seq 0 $LengthfileCdOneminusOne`
do

#il faut mettre les différentes options de commande sur une seule ligne sinon cela ne marche pas.
salmon quant -i $salmon_index \
--threads $nbThread \
--seqBias \
--gcBias \
-l A \
-1 ${arrayfastq[$i]}$extR1 \
-2 ${arrayfastq[$i]}$extR2 \
-o ${pathresults}/${arrayfastq[$i]} \
--validateMappings 

done


#désactivation de conda
conda deactivate

