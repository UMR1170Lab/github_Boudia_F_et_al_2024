# Documentation of the Single-Cell Analysis


##  From the fastq to a normalized matrix : the Single Cell Pipeline of Gustave Roussy  


Fastqs were handled using the Single Cell pipeline of Gustave Roussy (v1.3). All the documentation associated to this pipeline can be found here : https://github.com/gustaveroussy/single-cell. All the processing of the data from the quantification of reads to the generation of normalized counts matrix (and the selection of the number of PCA dimensions as well as the number of clusters for the generation of UMAP representations) were performed using this pipeline. For each sample, the pipeline has been launched and has allowed to select both a suitable number of PCA dimensions for umap representations and a reasonable number of clusters (the pipeline has already removed bad quality cells according to criteria provided in our paper). The table below gives the number of pca dimensions and the resolution of the clustering kept for each sample in the first part of our analysis.

 
| Sample | Number of PCA dimensions | Resolution for the clustering |
|--------|--------------------------|-------------------------------|
|CPN CTRL 1| 59 | 1.3 |
|----------|----|-----|
|CPN CTRL 2| 81 | 1.5 |
|----------|----|-----|
|CPN EG 1| 83 | 1.6 |
|--------|----|-----|
|CPN EG 2| 85 | 1.3 |
|--------|----|-----|
|NSG 441 |87  | 1.5 |
|--------|----|-----|
|NSG IL 514| 85 |1.3|
|----------|----|---|
|patient 1 (AML64) | 91 | 0.9 |
|------------------|----|-----|
|patient (GHER) | 95 | 0.9 |
|---------------|----|-----|
|patient (JER) | 75 | 0.1 |
|--------------|---- |----|
