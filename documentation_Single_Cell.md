# Documentation of the Single-Cell Analysis


##  From the fastq to a normalized matrix : the Single Cell Pipeline of Gustave Roussy  


Fastqs were handled using the Single Cell pipeline of Gustave Roussy (v1.3). All the documentation associated to this pipeline can be found here : https://github.com/gustaveroussy/single-cell. All the processing of the data from the quantification of reads to the generation of normalized counts matrix (and the selection of the number of PCA dimensions as well as the number of clusters for the generation of UMAP representations) were performed using this pipeline. For each sample, the pipeline has been launched and has allowed to select both a suitable number of PCA dimensions for umap representations and a reasonable number of clusters (the pipeline has already removed bad quality cells according to criteria provided in our paper). The table below gives the number of pca dimensions and the resolution of the clustering kept for each sample.  
