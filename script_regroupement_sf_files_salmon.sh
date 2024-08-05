#!/bin/bash


############################ Objectifs du script ###########################


## regrouper les fichiers de quantification .sf


#################### Variables ##################

#path du dossier contenant les fichiers de quantification .sf de Salmon.
path_quantif_file="/path/to/output/quantification/output_quantification_salmon" 

#path dossier où vont être regrouper les fichiers de quantification.

path_sf_files_groupir="/path/to/directory/regroupement_fichiers_sf_quantification_salmon"


################### code à exécuter #####################

#déplacements dans le dossier contenant les fichiers fastqs
cd $path_quantif_file


#tous les fichiers sf sont copiés dans un même dossier
for i in `ls`;do find $i -name '*.sf' -exec cp {} ${path_sf_files_groupir}/${i}.sf \; ;done




