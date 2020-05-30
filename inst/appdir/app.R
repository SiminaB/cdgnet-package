library(BiocManager)
options(repos = BiocManager::repositories())

library(CDGnet)

CDGnet:::.check_KEGG()
load(file.path(system.file("appdir", package="CDGnet"), "list_paths_KEGG.RData"))

data(DrugBank_targets)
data(drugs_PO_FDA_biomarkers)
data(drugs_PO_FDA_targets)
data(KEGG_cancer_paths_onc_long)
data(FDA_approved_drugs)
data(non_prop2prop)
data(prop2non_prop)
data(prop_non_prop)
data(Onc_df)

shiny::runApp(CDGnet::cdgnetApp())
