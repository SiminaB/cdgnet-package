require(KEGGREST)
require(KEGGgraph)
require(org.Hs.eg.db)
require(dplyr)

.check_and_load_KEGG <- function(ask_and_install=FALSE) {
  final_out_path <- file.path(system.file(package="CDGnet"), "kegg_data")
  filepath <- file.path(final_out_path, "list_paths_KEGG.RData")

  if (!file.exists(filepath)) {
    cat("[CDGnet] KEGG pathway data has not been downloaded yet. This needs one time to be done to properly use the CDGnet application.\n")

    if (ask_and_install) {
      cat("[CDGnet] Do you want to download KEGG pathway data now?")
      if (utils::menu(c("Yes", "No")) == 1) {
        .download_and_process_KEGG()
      } else {
        cat("[CDGnet] You chose not to download KEGG pathway data. The CDGnet app will not function properly. ")
        cat("You can download data at a later time using function 'CDGnet:::.download_and_process_KEGG()'\n")
        stop()
      }
    } else {
      cat("You can download data using function 'CDGnet:::.download_and_process_KEGG()'")
      stop()
    }
  }

  message("[CDGnet] Loading KEGG pathway data")
  load(filepath)
}

.download_and_process_KEGG <- function(basedir=tempdir()) {
  kgml_path <- file.path(basedir, "_KGML")
  if (!file.exists(kgml_path)) {
    dir.create(kgml_path)
  }

  csv_path <- file.path(basedir, "_csv")
  if (!file.exists(csv_path)) {
    dir.create(csv_path)
  }

  final_out_path <- file.path(system.file(package="CDGnet"), "kegg_data")
  if (!file.exists(final_out_path)) {
    dir.create(final_out_path)
  }

  .get_all_KEGG_hsa_pathways(kgml_path)
  .parse_KGML(in_path=kgml_path, out_path=csv_path)
  res <- .bind_KEGG_files(csv_path)
  list_paths_KEGG <- res$list_paths_KEGG
  save(list_paths_KEGG, file=file.path(final_out_path, "list_paths_KEGG.RData"))

  #KEGG_cancer_paths_onc_long <- .parse_KEGG_oncogene_info(connections_KEGG)
  #.preprocess_KEGG_objects(KEGG_cancer_paths_onc_long,
  #                         list_paths_KEGG,
  #                         final_out_path)
}

.preprocess_KEGG_objects <- function(KEGG_cancer_paths_onc_long,
                                     list_paths_KEGG_stand,
                                     out_path) {
  message("[cdgnet] Saving KEGG info")

  ##change a couple of the pathway names in KEGG to be more recognizable
  KEGG_cancer_paths_onc_long$Name <-
    as.character(KEGG_cancer_paths_onc_long$Name)
  KEGG_cancer_paths_onc_long <- KEGG_cancer_paths_onc_long %>%
    mutate(Name =
             replace(Name,
                     Name == "Glioma",
                     "Glioblastoma"))
  KEGG_cancer_paths_onc_long <- KEGG_cancer_paths_onc_long %>%
    mutate(Name =
             replace(Name,
                     Name == "Urothelial carcinoma",
                     "Bladder cancer"))
  ##for some reason, "Gastric cancer" did not get downloaded - take it out!
  KEGG_cancer_paths_onc_long <- KEGG_cancer_paths_onc_long %>%
    filter(Name != "Gastric cancer")

  ##take subset of KEGG pathways that only correspond to cancer pathways AND
  ##make the list of KEGG pathways to have as the keys the pathway names NOT the pathway IDs
  ##first make an id2name vector
  KEGG_unique_paths <- unique(KEGG_cancer_paths_onc_long[,1:2])
  KEGG_path_id2name <- as.character(KEGG_unique_paths$Name)
  names(KEGG_path_id2name) <- as.character(KEGG_unique_paths$KEGG_id)
  ##take subset of KEGG pathways that are ONLY cancer pathways
  list_paths_KEGG <- list_paths_KEGG[names(KEGG_path_id2name)]
  names(list_paths_KEGG) <- KEGG_path_id2name[names(list_paths_KEGG)]

  save(KEGG_cancer_paths_long_onc, file=file.path(out_path, "KEGG_cancer_paths_long_onc.RData"))
  save(list_paths_KEGG, file=file.path(out_path, "list_paths_KEGG.RData"))
}

.parse_KEGG_oncogene_info <- function(KEGG_cancer_paths_onc) {
  message("[cdgnet] Cleaning up KEGG info")

  KEGG_cancer_paths_onc <- KEGG_cancer_paths_onc[,1:3]
  KEGG_cancer_paths_onc$KEGG_id <- as.character(KEGG_cancer_paths_onc$KEGG_id)
  KEGG_cancer_paths_onc$Name <- as.character(KEGG_cancer_paths_onc$Name)
  KEGG_cancer_paths_onc$Oncogen <- as.character(KEGG_cancer_paths_onc$Oncogen)

  onc_list <- strsplit(as.character(KEGG_cancer_paths_onc$Oncogenes), ";")
  names(onc_list) <- as.character(KEGG_cancer_paths_onc$KEGG_id)

  KEGG_cancer_paths_onc_long <- data.frame(KEGG_id = character(),
                                           Name = character(),
                                           Oncogene = character())

  for (KEGG_id in names(onc_list)) {
    onc_list_id <- onc_list[[KEGG_id]]
    name_id <-
      KEGG_cancer_paths_onc$Name[KEGG_cancer_paths_onc$KEGG_id==KEGG_id]
    KEGG_cancer_paths_onc_long <- rbind(KEGG_cancer_paths_onc_long,
                                        data.frame(KEGG_id=KEGG_id,
                                                   Name=name_id,
                                                   Oncogene=onc_list_id))
  }
  KEGG_cancer_paths_onc_long
}

.bind_KEGG_files <- function(in_path) {
  message("[cdgnet] Processing KEGG graph")

  all_KEGG <- list.files(in_path)
  connections_KEGG <- data.frame()
  list_paths_KEGG <- list()
  for(KEGG_path in all_KEGG)
  {
    KEGG_path_df <- read.csv(paste(in_path,
                                   KEGG_path, sep="/"))
    connections_KEGG <- rbind(connections_KEGG,
                              KEGG_path_df)
    KEGG_path_name <- gsub(".csv","",KEGG_path)
    list_paths_KEGG[[KEGG_path_name]] <- KEGG_path_df[,-1]
  }

  ##get rid of first column
  connections_KEGG <- connections_KEGG[,-1]
  res <- list(connections_KEGG=connections_KEGG,
              list_paths_KEGG=list_paths_KEGG)

  res
}

.parse_KGML <- function(in_path, out_path) {
  message("[cdgnet] Parsing downloaded KGML files")

  all_kgml <- list.files(in_path)
  for (file_kgml in all_kgml) {
    file_path_kgml <- paste(in_path, file_kgml, sep="/")

    ##convert KGML to data frame
    current_df <- parseKGML2DataFrame(file_path_kgml)

    if(nrow(current_df) > 0)
    {
      current_df$from <- as.character(current_df$from)
      current_df$to <- as.character(current_df$to)
      ##get all KEGG IDs for nodes
      all_nodes <- unique(c(current_df[,1], current_df[,2]))

      ##first get rid of the "hsa:" part in front of the IDs
      all_nodes_gid <- gsub("hsa:","",all_nodes)
      names(all_nodes_gid) <- all_nodes
      if (sum(!is.na(all_nodes_gid)) > 0) {
        length(unique(all_nodes))
        length(unique(all_nodes_gid))
        head(all_nodes)
        head(all_nodes_gid)

        gid_to_name <- as.character(org.Hs.egSYMBOL)
        all_nodes_gene_symbol <- gid_to_name[all_nodes_gid]
        head(all_nodes_gene_symbol)
        length(all_nodes_gene_symbol)
        all_nodes_gene_symbol <- gid_to_name[all_nodes_gid]

        ##now put in gene names in the data frame
        ##only keep KEGG IDs that do get mapped
        ##(the other IDs are probaly things like pathways _not_ genes/proteins)
        current_df <- dplyr::filter(current_df,
                                    from %in% names(all_nodes_gid),
                                    to %in% names(all_nodes_gid))

        current_df$from <- all_nodes_gene_symbol[all_nodes_gid[current_df$from]]
        current_df$to <- all_nodes_gene_symbol[all_nodes_gid[current_df$to]]
        # head(current_df)
        #
        # graph <- graph_from_data_frame(current_df)
        # V(graph)$size <- 8
        # V(graph)$label.cex <- 0.6

        # png_name <- paste("../All_human_pathways_figs",gsub("(.xml)|(.kgml)", ".png",file_kgml),sep="/")
        #
        # png(filename=png_name)
        # try(plot(graph,
        #          edge.arrow.size=0,
        #          layout=layout_with_kk))
        # dev.off()

        csv_name <- paste(out_path, gsub("(.xml)|(.kgml)", ".csv", file_kgml), sep="/")
        write.csv(x=current_df, file=csv_name)

      }
    }
  }
}

.get_all_KEGG_hsa_pathways <- function(out_path) {
  all_path_KEGG <- keggList("pathway", organism="hsa")

  message("[cdgnet] Downloading pathway info from KEGG (this may take a few minutes)")

  suppressWarnings(for(path in names(all_path_KEGG)) {
    path_name <- gsub("path:","",path)
    path_kgml <- keggGet(path, "kgml")
    filename <- paste0(path_name, ".kgml")
    write(path_kgml, file.path(out_path, filename))
  })
}
