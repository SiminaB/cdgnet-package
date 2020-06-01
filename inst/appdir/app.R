library(BiocManager)
options(repos = BiocManager::repositories())

library(CDGnet)

if (!file.exists("list_paths_KEGG.RData")) {
  stop("KEGG data not found. Download with function CDGnet::download_and_process_KEGG()")
}

load("list_paths_KEGG.RData")

cat("list_paths_KEGG names:")
print(names(list_paths_KEGG))

data(DrugBank_targets)
data(drugs_PO_FDA_biomarkers)
data(drugs_PO_FDA_targets)
data(KEGG_cancer_paths_onc_long)
data(FDA_approved_drugs)
data(non_prop2prop)
data(prop2non_prop)
data(prop_non_prop)
data(Onc_df)

CDGnet::cdgnetApp()
