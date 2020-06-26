library(BiocManager)
options(repos = BiocManager::repositories())

library(CDGnet)

CDGnet:::.check_KEGG(package=FALSE)
kegg_path <- CDGnet:::.keggFile(package=FALSE)
cdgnetApp(kegg_path)
