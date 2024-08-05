#!/usr/bin/bash
#SBATCH --partition=mediumq
#SBATCH --time=24:00:00
#SBATCH --output=rapport_nfcore_pipeline_atacseq.out
#SBATCH --mem=90G

########################### chargement des modules nécessaires au lancement du pipeline nextflow ###########################



module load java/1.8.0_181
#chargement de singularity
module load singularity
#chargement de nextflow
module load nextflow/20.04.1 

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh


#path vers fichier csv contenant le design

pathdesign="/path/to/design.csv"

pathresults="/path/to/results_v1"

pathconfigfile="/path/to/config_file.config"

Djava_tmpdir="/path/to/java_tmp"

#-Djava.io.tmpdir=$Djava_tmpdir

########################################### lancement du pipeline nextflow

nextflow run /mnt/beegfs/pipelines/nf-core-atacseq/1.2.1/main.nf \
--input $pathdesign \
--genome hg38 \
--outdir $pathresults \
-resume \
-c $pathconfigfile \
--igenomes_base /path/to/hg38_genome_for_nfcore
