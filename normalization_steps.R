set.seed(seed = 1)
    
#unormalised_datasets corresponds to an R named list with all the datasets that should be integrated. 


#creation of the variable alltogether which correspond to the integration of the unnormalised datasets
alltogether <- unnormalised_datasets

alltogether <- lapply(X = alltogether, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = alltogether)

alltogether <- lapply(X = alltogether, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})


alltogether_anchors <- FindIntegrationAnchors(object.list = alltogether, anchor.features = features, reduction = "rpca")


# this command creates an 'integrated' data assay
alltogether.combined <- IntegrateData(anchorset = alltogether_anchors)

