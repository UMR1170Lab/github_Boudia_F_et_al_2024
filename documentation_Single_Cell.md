# Documentation of the Single-Cell Analysis


##  From the fastq to a normalized matrix : the Single Cell Pipeline of Gustave Roussy  


Fastqs were handled using the Single Cell pipeline of Gustave Roussy (v1.3). All the documentation associated to this pipeline can be found here : https://github.com/gustaveroussy/single-cell. All the processing of the data from the quantification of reads to the generation of normalized counts matrix (and the selection of the number of PCA dimensions as well as the number of clusters for the generation of UMAP representations) were performed using this pipeline. For each sample, the pipeline has been launched and has allowed to select both a suitable number of PCA dimensions for umap representations and a reasonable number of clusters (the pipeline has already removed bad quality cells according to criteria provided in our paper). The table below gives the number of pca dimensions and the resolution of the clustering kept for each sample in the first part of our analysis.

 
| Sample | Number of PCA dimensions | Resolution for the clustering |
|--------|--------------------------|-------------------------------|
|CPN CTRL 1| 59 | 1.3 |
|CPN CTRL 2| 81 | 1.5 |
|CPN EG 1| 83 | 1.6 |
|CPN EG 2| 85 | 1.3 |
|NSG 441 |87  | 1.5 |
|NSG IL 514| 85 |1.3|
|patient 1 (AML64) | 91 | 0.9 |
|patient (GHER) | 95 | 0.9 |
|patient CONECT| 73 |1.4|

An example of the files used to launch the pipeline is provided in this branch : run_CPN_NEG_disomique_IPS_CTRL_jour13_sans_stroma.sh and CPN_NEG_disomique_IPS_CTRL_jour13_sans_stroma.yaml

For our study, only cells carrying the ETO2-GLIS2 mutation were of interest in patient samples. For this reason, we have identified cells of the immune infiltrate in these samples and removed them (removal of these cells in the unormalized matrix without bad quality cells). The single-cell pipeline was launched again starting from the unormalized matrix without immune cells and bad quality cells. The table below gives the number of pca dimensions and the resolution of the clustering kept for each patient's sample without cells of the immune infiltrate.

| Sample | Number of PCA dimensions | Resolution of the clustering |
|--------|--------------------------|------------------------------|
|patient 1 (AML64) | 77 | 0.8 |
|patient 2 (GHER) | 91 | 0.8 |
|patient 3 (CONECT) | 93 | 0.8 |

## Integration of all the Single-Cell samples

Once all samples that we wanted to integrate have been properly handled (generation of a normalized matrix and removal of cells belonging to the immune infiltrate in patient's samples), we used the RPCA integration method provided in the R package Seurat to integrate them. The script that has been used to perform the integration is normalization_steps.R. For UMAp representations, we have decided to keep 20 dimensions of PCA. The clustering resolution kept was of 0.5.







 
