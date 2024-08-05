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
|patient C  (AML64) | 91 | 0.9 |
|patient A (GHER) | 95 | 0.9 |
|patient B  CONECT| 73 |1.4|

An example of the files used to launch the pipeline is provided in this branch : run_CPN_NEG_disomique_IPS_CTRL_jour13_sans_stroma.sh and CPN_NEG_disomique_IPS_CTRL_jour13_sans_stroma.yaml

For our study, only cells carrying the ETO2-GLIS2 mutation were of interest in patient samples. For this reason, we have identified cells of the immune infiltrate in these samples (annotations generated with the R package SingleR based on the dataset of Novershtern N et al. (2011)) and removed them (removal of these cells in the unormalized matrix without bad quality cells). The single-cell pipeline was launched again starting from the unormalized matrix without immune cells and bad quality cells. The table below gives the number of pca dimensions and the resolution of the clustering kept for each patient's sample without cells of the immune infiltrate. 

| Sample | Number of PCA dimensions | Resolution of the clustering |
|--------|--------------------------|------------------------------|
|patient C (AML64) | 77 | 0.8 |
|patient A (GHER) | 91 | 0.8 |
|patient B (CONECT) | 93 | 0.8 |

## Integration of all the Single-Cell samples

Once all samples that we wanted to integrate have been properly handled (generation of a normalized matrix and removal of cells belonging to the immune infiltrate in patient's samples), we used the RPCA integration method provided in the R package Seurat to integrate them. The script that has been used to perform the integration is normalization_steps.R.


## testing specific signature : use of the R function projection_signature_integration_v2

This function allows to sort cells in two category : cells in which we can detect genes of the signature in a significant higher proportion than expected randomly and cells in whih the proportion of detected genes of the signature is not higher than expected randomly. To use this function, the following code should be previously run :

```
require(purrr)

require(Seurat)

source('/path/in/your/computer/of/pvaldistri_v3.R')

source('/path/in/your/computer/of/association_number_condition.R')

source('/path/in/your/computer/of/df_rep_of_each_level_a_factor_in_each_level_another_factor.R')

source('/path/in/your/computer/of/IndiceRepSignaturePval_v5.R')


source('/path/in/your/computer/of/projection_signature_integration_v2.R')

```

For example to test the EG signature as done in the article, run :

```

gene_up_EG_Seattle <- read.table("/path/in/your/computer/of/the/file/geneup_in_EG_compare_to_others_LFCsup1.gmx",h=T)

EGLFCupsup1 <- gene_up_EG_Seattle$EGLFCupsup1[-1]


projection_signature_integration_v2(selected_integration,
                                    normalized_individual_datasets,
                                    interest_biological_signature = EGLFCupsup1[1:500],
                                    path_signature_statut_cells_save_rds = "/path/in/your/computer/of/the/directory/where/you/want/to/store/your/output/file_that_will_store_the_classification_for_each_cells.rds" ,
                                    name_biological_signature = "EG_LFCsup1_top500", 
                                    minimal_prop_cells = 0.02,
                                    number_random_signa = 100,
                                    threshold_pval = 0.01,
                                    used_random_seed = 100,
                                    path_pdf_figures = "/path/in/your/computer/of/the/directory/where/you/want/to/store/your/output/file_containing_the_results_representation.pdf",
                                    width_pdf = 30,
                                    height_pdf = 30)


```


 



 

