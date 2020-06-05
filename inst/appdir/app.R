library(BiocManager)
options(repos = BiocManager::repositories())

library(CDGnet)

list_paths_KEGG_file <- file.path(getwd(), "list_paths_KEGG.RData")

if (!file.exists(list_paths_KEGG_file)) {
  stop("KEGG data not found. Download with function CDGnet::download_and_process_KEGG()")
}

CDGnet::cdgnetApp(list_paths_KEGG_file)
