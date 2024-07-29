#!/bin/bash

########################################################################
## Single-cell RNA-seq script to launch single-cell RNA-seq pipeline
##
## using: sbatch /mnt/beegfs/userdata/m_aglave/pipeline/run.sh
##
########################################################################

## JOB PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#SBATCH --job-name=pipeline_sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=mediumq
#SBATCH --output=individual_analysis_CPN_NEG_disomique_IPS_CTRL_jour13_sans_stroma_Clust_Markers_Annot_GE_Cerebro.out

source /path/to/conda.sh
conda activate /path/to/.environnement_conda/conda_pipelinesingle_cell_marine
module load singularity

#print environment tools versions
python --version
snakemake --version
singularity --version

#parameters
path_to_configfile="/path/to/CPN_NEG_disomique_IPS_CTRL_jour13_sans_stroma.yaml"
path_to_pipeline="/path/to/pipelines/single-cell/1.3"

#launch
snakemake --profile ${path_to_pipeline}/profiles/slurm \
-s ${path_to_pipeline}/Snakefile \
--configfile ${path_to_configfile}


conda deactivate
