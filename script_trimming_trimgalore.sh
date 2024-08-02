#!/bin/bash
#SBATCH --job-name=trimm_galore               # Nom du Job
#SBATCH --ntasks=1                       # Nombre de Tasks : 1
#SBATCH --cpus-per-task=12                # Allocation de 12 CPU par Task
#SBATCH --time=24:00:00
#SBATCH --mem=35G
#SBATCH --partition=mediumq
#SBATCH --output=rapport_trimming.out 


#Date d'écriture du script : 26 juin 2021

##################### fonction du script ###############

##virer les séquences d'adaptateurs des fastqs de Loelia concernant NPM-ALK


################ Variables utilisés dans le script ###############

## path pour activer conda
pathconda="/path/to/conda.sh"

##Environnment conda qui va être utilisé pour ce script.
condatrimmgalore="envi_remove_adaptator"

### Dossier qui contient les fastqs non trimmés
fastqfile_before_trimming="/path/to/symbolic_link_fastq"


### Variable qui permet d'indiquer si on a affaire à du single-end ou à du paired-end
SE="PE"

## variable qui donne le suffixe des fichiers R1
suffixeR1="_R1_001.fastq.gz"

#Variable qui donne le path dans lequel vont se trouver les fichiers trimm_galore.
pathrepertorytrimmgalore="/path/to/directory/fastq_after_trimming"


NbThread=11

######################## lancement du script ########################


#activation de conda
source $pathconda

#activation de l environnement conda contenant trimm galore

conda activate $condatrimmgalore


#Création du dossier qui va contenir les fastqs une fois que le trimming sera terminé.

mkdir ${pathrepertorytrimmgalore}/output_fastq_trimmgalore


# Déplacement dans le dossier contenant les fichiers fastqs non trimmés

cd $fastqfile_before_trimming


#Si les fichiers fastqs sont single-end
if [ $SE = "SE" ]; then

echo "trimming pour du single-end"

while read FASTQ; do

fastqpref=$(echo $FASTQ | sed "s/$suffixeR1//")


#préfixe du fichier sans le chemin d'accès qui lui est associé
basefastqpref=$(basename $fastqpref)


echo $fastqpref
echo $basefastqpref


echo -e "Beginning of the removal of the adapters sequences and of poor qualities bases"

trim_galore --output_dir ${pathrepertorytrimmgalore}/output_fastq_trimmgalore \
--basename $basefastqpref \
--phred33 --fastqc --quality 20 --cores NbThread $fastqpref$suffixeR1 > ${pathrepertorytrimmgalore}/output_fastq_trimmgalore/${basefastqpref}.log 2>&1

done < <(
find ${fastqfile_before_trimming}/ -name '*R1*.fastq.gz'
);

elif [ $SE = "PE" ];then

echo "trimming pour du paired-end"

#définition de la valeur de la variable qui donne le sufix des fichiers R2
suffixeR2=$(echo $suffixeR1 | sed "s/R1/R2/")

while read FASTQ; do

fastqpref=$(echo $FASTQ | sed "s/$suffixeR1//")


#préfixe du fichier sans le chemin d'accès qui lui est associé
basefastqpref=$(basename $fastqpref)


echo $fastqpref$suffixeR1 
echo $fastqpref$suffixeR2
echo $basefastqpref


echo -e "Beginning of the removal of the adapters sequences and of poor qualities bases"

trim_galore --output_dir ${pathrepertorytrimmgalore}/output_fastq_trimmgalore \
--basename $basefastqpref \
--phred33 \
--fastqc \
--quality 20 \
--paired \
--cores $NbThread \
$fastqpref$suffixeR1 $fastqpref$suffixeR2 > ${pathrepertorytrimmgalore}/output_fastq_trimmgalore/${basefastqpref}.log 2>&1

done < <(
find ${fastqfile_before_trimming}/ -name '*R1*.fastq.gz'
);

else

echo "je ne sais pas si les fichiers fastqs sont stranded ou unstranded"

fi

#désactivation de conda 
conda deactivate


conda activate multiqcenv

# déplacement dans le répertoire contenant les fastqs trimmés.
cd $pathrepertorytrimmgalore

multiqc .

conda deactivate


