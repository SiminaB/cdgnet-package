##function that maps to various attributes
# source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")
##see also http://michael.hahsler.net/SMU/LearnROnYourOwn/code/igraph.html for examples

#' function that inputs a list of genes and outputs a vector matching the names to the standardized names
#' for the genes that were matched up using DGI (and just returns the original names for the ones that were not)
#'
#' @param genes list of gene names to standardize
#' @return list of gene names
#' @import rDGIdb
#'
gene_names_standard <- function(genes)
{
  stg <- searchTermSummary(queryDGIdb(genes))
  ## make list with standardized gene names as names, initial gene names as entries
  gene_names_list <- as.list(as.character(stg$SearchTerm))
  names(gene_names_list) <- as.character(stg$Matches)
  ## now split on comma
  gene_names_sep_list <- lapply(gene_names_list,
                                function(g){strsplit(g,", ")[[1]]})
  ## now create vector for id matching
  ## first repeat the standard names the correct number of times
  stand_names <- rep(names(gene_names_sep_list), sapply(gene_names_sep_list, length))
  name2stand_name <- stand_names
  names(name2stand_name) <- unlist(gene_names_sep_list)

  ## now print out all the genes that do not get mapped to standard names
  print("Genes that do not get mapped to standard names:")
  not_mapped <- setdiff(genes, names(name2stand_name))
  print(not_mapped)
  ## remove NA
  not_mapped <- not_mapped[!is.na(not_mapped)]

  ## for all these genes, add in the IDs at the end (except for NAs)
  name2stand_name <- c(name2stand_name, not_mapped)
  names(name2stand_name) <- c(unlist(gene_names_sep_list), not_mapped)

  name2stand_name
}

#' function that adds other genes to the molecular profiling (MP) data frame
#'
#' @param MP_df data.frame containing molecular profile
#' @param add_genes list of genes to add to profile
#' @return data frame
#'
MP_add_genes <- function(MP_df, add_genes)
{
  keep_genes_not_MP <- setdiff(add_genes, MP_df$Gene_protein)
  MP_plus <- MP_df
  if(length(keep_genes_not_MP) > 0)
  {
    MP_plus <- rbind(MP_df, data.frame(Gene_protein = keep_genes_not_MP,
                                       Data_type = NA,
                                       Alteration = NA))
  }
  MP_plus
}

#' function to remove wild type alterations from an MP data frame
#'
#' @param MP_df data frame with molecular alterations
#' @param wild_type ways in which the wild type can be coded in MP_df
#' @return data frame
#'
remove_WT <- function(MP_df,
                      wild_type = c("wild type","wildtype","wild","wt"))
{
  MP_df[!(tolower(MP_df$Alteration) %in% wild_type),]
}

#' function that checks whether any of the drugs have biomarkers that are wild-type
#' if any of the biomarkers are wild-type, then check the MP report
#' if they are not wild-type there, then remove them from the list
#'
#' @param drugs_bio data frame with drugs and biomarkers
#' @param MP_df data frame with molecular alterations
#' @return data frame
#'
#' @importFrom dplyr filter
#' @import magrittr
#'
check_WT_biomarkers <- function(drugs_bio,MP_df){

  ##if any of these are biomarkers for "wild type", check against that in MP_df data
  biomarkers_WT <- unique(drugs_bio$Gene.Protein[drugs_bio$Alteration == "wild type"])

  if(length(biomarkers_WT)>=1)
  {
    ##what is the alteration for these biomarkers in the MP_df report? If it is not "wild type" or "WT," remove it from the list
    biomarkers_WT_MP_alt <- MP_df %>% filter(Gene_protein %in% biomarkers_WT)
    biomarkers_WT_MP_alt <-
      biomarkers_WT_MP_alt[!(tolower(biomarkers_WT_MP_alt$Alteration) %in%
                               c("wild type","wildtype","wild","wt")),]
    keep_MP_bio <- setdiff(MP_df$Gene_protein, biomarkers_WT_MP_alt$Gene_protein)

    drugs_bio <- filter(drugs_bio, Gene.Protein %in% keep_MP_bio)
  }

  drugs_bio
}

#' function that calculates the shortest paths from inputs to everything else that
#' is in a pathway but not in MP
#' calls on calculate_shortest_path
calculate_shortest_path_from_inputs <- function(cancer_path,
                                                MP,
                                                inputs_in_path)
{
  ## convert the connections to an adjacency matrix
  g <- graph.data.frame(cancer_path[,c("from","to")])

  ##get the length of the shortest path from the inputs to all the other genes
  all_genes <- unique(c(as.character(cancer_path$from),
                        as.character(cancer_path$to)))

  ##only do the rest of this if at least one input is in all_genes - otherwise return NULL
  shortest_paths_from_inputs <- NULL

  if(sum(inputs_in_path %in% all_genes)!=0)
  {
    inputs <- unique(as.character(MP$Gene_protein))
    all_genes_no_inputs <- setdiff(all_genes, inputs)

    shortest_paths_from_inputs <- calculate_shortest_path(inputs_in_path,
                                                          all_genes_no_inputs,
                                                          g)
  }

  shortest_paths_from_inputs
}

#' function that calculates the shortest paths from some inputs to some outputs
#'
#' @param input character vectors
#' @param output character vector
#' @param igraph_obj is created by converting a data frame of edges using graph.data.frame
#' @return a matrix with the inputs as rows, the outputs as columns, and the shortest path length as entries
#'
calculate_shortest_path <- function(inputs, outputs, igraph_obj)
{
  shortest_paths_from_inputs <-
    matrix(nrow=length(inputs), ncol=length(outputs))
  rownames(shortest_paths_from_inputs) <- inputs
  colnames(shortest_paths_from_inputs) <- outputs

  for(input in inputs)
  {
    shortest_paths_from_inputs[input,] <-
      sapply(outputs,
             function(to,from,igraph_obj){length(attributes(shortest_paths(igraph_obj,from,to)$vpath[[1]])$names)},
             input,igraph_obj)
  }
  shortest_paths_from_inputs
}

#' function that pastes the shortest paths between some inputs and some outputs
#' (inputs_in_path are downstream of the oncogenes)
#'
#' @param input character vectors
#' @param output character vector
#' @param igraph_obj is created by converting a data frame of edges using graph.data.frame
#' @return a character vector with the paths marked with "->", separated by ";" if there are multiple paths
#'
paste_shortest_path <- function(inputs, inputs_in_path, outputs, igraph_obj)
{
  shorts <- rep(NA, length(outputs))
  names(shorts) <- outputs
  for(gene in outputs)
  {
    if(gene %in% inputs){
      shorts[gene] <- gene
    } else {
      gene_path_from_input <- list()
      for(input in inputs_in_path)
      {
        gene_path_from_input[[input]] <-
          paste(attributes(shortest_paths(igraph_obj, input, gene)$vpath[[1]])$names,
                sep="", collapse=" -> ")
      }
      gene_path_from_input_vect <- unlist(gene_path_from_input)
      gene_path_from_input_vect <-
        gene_path_from_input_vect[gene_path_from_input_vect != ""]
      shorts[gene] <- paste(gene_path_from_input_vect,
                            sep="",collapse=";")

    }
  }
  shorts
}

#' function that prettifies output of drugs_mat
#'
#' @param drug_mat a drug matrix
#' @return data frame
#'
pretty_drugs_mat <- function(drugs_mat)
{
  colnames(drugs_mat)[4:6] <- c("Predicted effect",
                                "Data type",
                                "Alteration")

  ##remove/simplify some things further
  drugs_mat <- drugs_mat[,c("Drug","Gene.Protein","Data type",
                            "Alteration","Path","Disease","Predicted effect")]
  drugs_mat
}

#' function to get nodes downstream from given nodes
#' starts from the entire pathway - calls calculate_shortest_path_from_inputs and
#' get_downstream_nodes
#'
get_downstream_from_inputs <- function(cancer_path,
                                       MP,
                                       inputs_in_path)
{
  ##get the length of the shortest path between the genes/proteins in
  ##inputs_in_path and all the other genes (except for the ones that are WT in the
  ##actual MP)
  shortest_paths_from_inputs <-
    calculate_shortest_path_from_inputs(cancer_path,
                                        MP,
                                        inputs_in_path)

  ##if this is not null, get downstream nodes, otherwise just return an empty vector
  downstream_nodes <- c()
  if(!is.null(shortest_paths_from_inputs))
  {
    ##keep genes that have a shortest path > 0 from any of the oncogenes that are also positive in the MP report
    ##get just these genes
    ##want genes downstream of _all_ the inputs
    ##first get list of all downstream things
    downstream_nodes <- get_downstream_nodes(inputs_in_path,
                                             shortest_paths_from_inputs)
  }
  downstream_nodes
}

#' function to get nodes downstream from given nodes (which are the inputs)
#'
#' @param inputs is a character vector
#' @param shortest_paths_from_inputs is matrix generated by calculate_shortest_path
#' @return output is a list of downstream nodes for each input node
#'
get_downstream_nodes <- function(inputs, shortest_paths_from_inputs)
{
  downstream_from_inputs <- list()
  for(input in inputs){
    downstream_from_inputs[[as.character(input)]] <-
      names(which(shortest_paths_from_inputs[input,]>0))
  }
  ##since want the shortest path from _any_ of the inputs, take the column sums
  sum_shortest_path <- colSums(shortest_paths_from_inputs)
  ##the ones that are downstream from any of the inputs have the sum > 0
  downstream_from_inputs <- names(which(sum_shortest_path>0))
  downstream_from_inputs
}

#' filter data frames based on the gene/protein being in keep_genes
#' keep only the unique drug names
#'
#' @param drugs_genes_df data.frame contains gene/protein names
#' @param keep_genes list of genes to filter or keep
#' @return data.frame
#'
#' @importFrom dplyr filter
#'
filter_drugs_select_genes <- function(drugs_genes_df,
                                      keep_genes)
{
  filter(drugs_genes_df,
                Gene.Protein %in% keep_genes)
}

#' function that creates the data frame specifying the edges from a data frame like drugs_PO_FDA_keep
#'
#' @param df_targeted data.drame
#' @param subtype default to sensitivity
#' @return data.frame
#'
make_df_edges_targeted <- function(df_targeted, subtype = "sensitivity")
{
  if(nrow(df_targeted)>0)
  {
    df_edges_targeted <- data.frame(from=as.character(df_targeted$Drug),
                                    to=as.character(df_targeted$Gene.Protein),
                                    subtype=subtype,
                                    direction=-1)
  } else {
    df_edges_targeted <- data.frame(from=character(),to=character(),
                                    subtype=character(),direction=character())
  }
  df_edges_targeted
}

#' function that combines all drugs from all categories
#'
#' @param Cat1 categeory 1 data.frame
#' @param Cat2, categeory 2 data.frame
#' @param Cat3, categeory 3 data.frame
#' @param Cat4 categeory 4 data.frame, can be null
#' @return merged data frame
#'
combine_drugs <- function(Cat1, Cat2, Cat3, Cat4 = NULL)
{
  drugs <- c()
  if(ncol(Cat1)>1)
  {
    drugs <- c(drugs, as.character(Cat1$Drug))
  }
  if(ncol(Cat2)>1)
  {
    drugs <- c(drugs, as.character(Cat2$Drug))
  }
  if(ncol(Cat3)>1)
  {
    drugs <- c(drugs, as.character(Cat3$Drug))
  }
  if(!is.null(Cat4))
  {
    if(ncol(Cat4)>1)
    {
      drugs <- c(drugs, as.character(Cat4$Drug))
    }
  }
  unique(drugs)
}

#' function that generates pathway plots for categories 3 and 4
#' @param edge is the matrix of edges
#' @param drugs character vectors
#' @param inputs character vectors
#' @return pathway plot
graph_pathways_cats_3_4 <- function(edges, drugs, inputs)
{
  graph <- graph_from_data_frame(edges)
  V(graph)$size <- 3
  V(graph)$label.cex <- 0.4
  V(graph)$color <- "#1f77b4"

  V(graph)$color[V(graph)$name %in% drugs] <- "#ff7f0e"
  V(graph)$color[V(graph)$name %in% inputs] <- "#9467bd"

  # ##get smallest distance between each vertex and each drug
  # ##only need to consider drugs that are vertex names
  # drugs <- intersect(drugs, V(graph)$name)
  # small_dist_drug <- calculate_shortest_path(drugs, V(graph)$name, graph)
  # ##replace 0 with number of vertices (probably means that there's no direct path)
  # small_dist_drug[small_dist_drug == 0] <- length(V(graph)$name)
  # ##get smallest distance to any drug
  # small_dist_any_drug <- apply(small_dist_drug, 2, min)
  # ##put in a very small number for the actual drugs
  # small_dist_any_drug[drugs] <- 1
  # write.csv(small_dist_drug, file="small_dist_drug.csv")
  # write.csv(small_dist_any_drug, file="small_dist_any_drug.csv")

  plot(graph,
       edge.arrow.size=0,
       layout=layout_with_fr,
       edge.arrow.size=0.5,
       vertex.label.cex=0.7,##0.5*map(1/small_dist_any_drug[V(graph)$name],c(1,1.5)),
       vertex.label.family="Helvetica",
       vertex.label.font=2,
       vertex.shape="circle",
       vertex.label.dist=1,
       vertex.size=4,##map(1/small_dist_any_drug[V(graph)$name],c(1,10)),
       vertex.label.color="black",
       edge.width=0.5)
}

#' function that makes data frame of nodes
#' @param assoc_df is the data frame of drug-gene/protein assocations
#' @return data frame
#'
make_node_df <- function(assoc_df)
{
  node_names <- unique(c(as.character(assoc_df$from),
                         as.character(assoc_df$to)))
  nodes_df <- data.frame(Name = node_names,
                         Group = "Gene/Protein")
  nodes_df$Group <- as.character(nodes_df$Group)
  nodes_df$Group[grep("\\(",as.character(nodes_df$Name))] <- "Drug"

  nodes_df
}

#' function that gets connections from one or more pathways which are within a list like the one of KEGG pathways
#' @param list_of_paths is the list of pathways
#' @param IDs is a vector of IDs for the names of the pathways that should be pulled out
#' @param subtype may specify "activation" or just be NULL, in which case everything is included
#' @return output is a data frame of connections
#'
#' @importFrom dplyr filter
#'
get_connections_paths <- function(list_of_paths, IDs, subtype=NULL)
{
  cancer_path <- c()

  for(ID in IDs)
  {
    cancer_path <- rbind(cancer_path,
                         list_of_paths[[ID]])
  }
  ##filter based on subtype (if applicable)
  if(!is.null(subtype))
  {
    cancer_path <- filter(cancer_path,
                                 subtype == subtype)
  }

  ##remove any duplicates
  cancer_path <- cancer_path[!duplicated(cancer_path),]

  cancer_path
}

#' get edges for categories 3 and 4
#'
#' @param Cat3 reactive objects corresponding to categories 3
#' @param Cat4 reactive objects corresponding to categories 4, Cat4 may be null
#' @return edges
#'
get_edges_cats_3_4 <- function(Cat3, Cat4 = NULL)
{
  all_edges <- as.data.frame(Cat3$edges)

  if(!is.null(Cat4))
  {
    all_edges <- rbind(all_edges,
                       as.data.frame(Cat4$edges))
  }

  drugs_mat <- Cat3$drugs_mat

  if(!is.null(Cat4))
  {
    drugs_mat <- rbind(drugs_mat,
                       as.data.frame(Cat4$drugs_mat))
  }

  ##drugs that are not used should be the set difference between everything in an edge and
  ##everything in drugs_mat that is either a drug or a gene/protein
  # drugs_excluded <- unique(setdiff(c(as.character(all_edges[,1]),
  #                                    as.character(all_edges[,2])),
  #                                  c(as.character(drugs_mat[,"Drug"]),
  #                                    as.character(drugs_mat[,"Gene or Protein"]))))
  #
  # all_edges <- all_edges[!(all_edges$Source %in% drugs_excluded) &
  #                          !(all_edges$Target %in% drugs_excluded),]
  #
  all_edges <- all_edges[!duplicated(all_edges),]

  all_edges
}

#' for category 3 and 4 recommendations, take out the drugs already recommended above
#' (i.e. recommended in categories 1 and 2 for category 3 and in categories 1, 2, or 3
#' for category 4)
#' if Type3_df == NULL, then this means just take out category 1 and 2 drugs + anything that has
#' (otherwise, also take out category 3 drugs)
#'
#' @param drugs_biomarkers_targets is the initial input data frame, from which known recommendations will be removed
#' @param Type1_df category 1 drug data frame
#' @param Type2_df category 2 drug data frame
#' @param Type3_df category 3 drug data frame
#' @return filtered data frame
#'
remove_known_recs <- function(drugs_biomarkers_targets, Type1_df, Type2_df, Type3_df = NULL)
{
  known_recs <- NULL

  if(ncol(Type1_df) > 1)
  {
    known_recs <- Type1_df[,c("Gene or Protein","Type","Alteration")]
  }
  if(ncol(Type2_df) > 1)
  {
    known_recs <- Type2_df[,c("Gene or Protein","Type","Alteration")]
  }
  if(ncol(Type1_df) > 1 & ncol(Type2_df) > 1)
  {
    known_recs <- rbind(Type1_df[,c("Gene or Protein","Type","Alteration")],
                        Type2_df[,c("Gene or Protein","Type","Alteration")])
  }
  if(!is.null(Type3_df))
  {
    if(ncol(Type3_df)>1)
    {
      known_recs <- rbind(known_recs[,c("Gene or Protein","Type","Alteration")],
                          Type3_df[,c("Gene or Protein","Type","Alteration")])
    }
  }
  if(!is.null(known_recs))
  {
    drugs_biomarkers_targets <-
      drugs_biomarkers_targets[!(drugs_biomarkers_targets$Gene.Protein %in% known_recs[,1]),]
  }
  drugs_biomarkers_targets
}

#' function to add in the shortest path to the drugs_biomarkers_targets data frame
#' @importFrom dplyr mutate
#'
add_shortest_path <- function(drugs_biomarkers_targets,
                              inputs, inputs_in_path,
                              cancer_path)
{
  ##convert the connections to an adjacency matrix
  g <- graph.data.frame(cancer_path[,c("from","to")])

  ##get shortest path from input genes to keep_genes
  shorts <- paste_shortest_path(inputs, inputs_in_path,
                                unique(drugs_biomarkers_targets$Gene.Protein),
                                g)

  mutate(drugs_biomarkers_targets, Path=shorts[Gene.Protein])
}


#' Parse the result category into the ui json object
#'
#' @param input molecular profile
#' @param fda fda approved drugs (cat 1)
#' @param fda_other fda approved drugs for other cancer types ( cat 2)
#' @param rec drug recommendations
#' @return json object
#'
dataParser <- function(input, fda, fda_other, rec) {

  input <- distinct(input)
  fda <- distinct(fda)
  fda_other <- distinct(fda_other)
  rec <- distinct(rec)

  ###
  ## groups denotes the id and title for each section in the sankey diagram.
  ###
  json <- list(
    "links"=list(),
    "nodes"=list(),
    "groups"=list(
      "mp"=list(
        "type"="process",
        "title"="Molecular Profile",
        "bundle"=NULL,
        "id"="mp",
        "nodes"=c(),
        "def_pos"=NULL
      ),
      "fda"=list(
        "type"="process",
        "title"="FDA Approved Drugs",
        "bundle"=NULL,
        "id"="fda",
        "nodes"=c(),
        "def_pos"=NULL
      ),
      "it"=list(
        "type"="process",
        "title"="Inferred Targets",
        "bundle"=NULL,
        "id"="it",
        "nodes"=c(),
        "def_pos"=NULL
      ),
      "rec"=list(
        "type"="process",
        "title"="Recommended Therapies",
        "bundle"=NULL,
        "id"="rec",
        "nodes"=c(),
        "def_pos"=NULL
      )
    ),
    "order"=list(),
    "alignLinkTypes"=FALSE
  )

  link <- list(
    "color"="rgb(204, 235, 197)",
    "source"=NULL,
    "value"=10,
    "type"=NULL,
    "target"=NULL,
    "opacity"=1,
    "time"="*",
    "title"=NULL,
    "id"=NULL
  )

  node <- list(
    "bundle"=NULL,
    "title"=NULL,
    "visibility"="visible",
    "def_pos"=NULL,
    "id"=NULL,
    "style"="process",
    "direction"="r"
  )

  node_count <- 1
  link_count <- 1

  mp_nodes <- c()
  fda_nodes <- c()
  it_nodes <- c()
  rec_nodes <- c()

  input2 <- aggregate(input, by=list(input$Gene_protein), FUN=paste)

  for (i in rownames(input2)) {
    row <- input2[i, ]
    mp_node <- node
    mp_node$title <- row$Gene_protein
    mp_node$id <- paste0("mp^", row$Gene_protein)
    mp_nodes <- c(mp_nodes, mp_node$id)

    json$nodes[[node_count]] <- mp_node
    node_count <- node_count + 1

    json$groups$mp$nodes <- c(json$groups$mp$nodes, paste0("mp^", row$Gene_protein))
  }

  if(nrow(fda) >= 1 && !length(fda$Note) >= 1 ) {

    unique_fda <- aggregate(fda, by=list(fda$Drug, fda$`Gene or Protein`), FUN=paste)


    for (i in rownames(unique_fda)) {
      row <- unique_fda[i, ]

      fda_node <- node
      fda_node$title <- row$Drug
      fda_node$id <- paste0("fda^", row$Drug)
      fda_nodes <- c(fda_nodes, fda_node$id)

      json$nodes[[node_count]] <- fda_node
      node_count <- node_count + 1

      json$groups$fda$nodes <- c(json$groups$fda$nodes, paste0("fda^", row$Drug))

      fda_link <- link
      fda_link$source <- paste0("mp^", row[["Gene or Protein"]])
      fda_link$target <- paste0("fda^", row$Drug)
      fda_link$type <- "FDA approved Drugs"
      fda_link$title <- paste0("alteration: ", row$Alteration)
      fda_link$id <- paste0(fda_link$source, " -> ", fda_link$target)

      json$links[[link_count]] <- fda_link
      link_count <- link_count + 1
    }
  }

  if(nrow(fda_other) >= 1 && !length(fda_other$Note) >= 1) {
    print(fda_other)
    unique_fda_other <- aggregate(fda_other, by=list(fda_other$Drug, fda_other$`Gene or Protein`), FUN=paste)


    for (i in rownames(unique_fda_other)) {
      row <- unique_fda_other[i, ]

      fda_other_node <- node
      fda_other_node$title <- row[["Group.1"]]
      fda_other_node$id <- paste0("fda^", row[["Group.1"]], " (other tumor type)")
      fda_nodes <- c(fda_nodes, fda_other_node$id)

      json$nodes[[node_count]] <- fda_other_node
      node_count <- node_count + 1

      json$groups$fda$nodes <- c(json$groups$fda$nodes,
                                 paste0("fda^", row[["Group.1"]], " (other tumor type)"))

      fda_link <- link
      fda_link$source <- paste0("mp^", row[["Group.2"]])
      fda_link$target <- paste0("fda^", row[["Group.1"]], " (other tumor type)")
      fda_link$type <- "FDA approved Drugs"
      fda_link$title <- paste0("alteration: ", row$Alteration, "\n <br> tumor drug is approved for: ", row[["Tumor in which it is approved"]])
      fda_link$id <- paste0(fda_link$source, " -> ", fda_link$target)

      json$links[[link_count]] <- fda_link
      link_count <- link_count + 1
    }
  }

  if(nrow(rec) >= 1) {
    unique_rec <- aggregate(rec, by=list(rec$Drug, rec$`Gene or Protein`), FUN=paste)

    for (i in rownames(unique_rec)) {
      row <- unique_rec[i, ]

      # rec_node <- node
      # rec_node$title <- row$Drug
      # rec_node$id <- paste0("rec^", row$Drug)

      # json$nodes[[node_count]] <- rec_node
      # node_count <- node_count + 1

      # json$groups$rec$nodes <- c(json$groups$rec$nodes, paste0("rec^", row$Drug))

      fda_link <- link
      fda_link$source <- paste0("it^", row[["Group.2"]])
      fda_link$target <- paste0("rec^", row[["Group.1"]])
      fda_link$type <- "Recommended therapies"
      fda_link$title <- paste0("alteration: ", row$Alteration, "\n <br> tumor drug is approved for: ", row[["Tumor in which it is approved"]],
                               "\n <br> pathway: ", row$Path)
      fda_link$id <- paste0(fda_link$source, " -> ", fda_link$target)

      json$links[[link_count]] <- fda_link
      link_count <- link_count + 1
    }

    for (i in unique(rec[["Drug"]])) {
      rec_node <- node
      rec_node$title <- i
      rec_node$id <- paste0("rec^", i)
      rec_nodes <- c(rec_nodes, rec_node$id)

      json$nodes[[node_count]] <- rec_node
      node_count <- node_count + 1

      json$groups$rec$nodes <- c(json$groups$rec$nodes, paste0("rec^", i))
    }

    for (i in unique(rec[["Gene or Protein"]])) {
      rec_node <- node
      rec_node$title <- i
      rec_node$id <- paste0("it^", i)
      it_nodes <- c(it_nodes, rec_node$id)

      json$nodes[[node_count]] <- rec_node
      node_count <- node_count + 1

      json$groups$it$nodes <- c(json$groups$it$nodes, paste0("it^", i))
    }

    for (i in unique(rec[["Path"]])) {

      paths <- strsplit(i, "->")

      rec_link <- link
      rec_link$source <- paste0("mp^", trimws(paths[[1]][1]))
      rec_link$target <- paste0("it^", trimws(paths[[1]][length(paths[[1]])]))
      rec_link$type <- "Pathway"
      rec_link$title <- i
      rec_link$id <- i

      json$links[[link_count]] <- rec_link
      link_count <- link_count + 1
    }
  }

  gmp <- json$groups$mp
  gfda <- json$groups$fda
  git <- json$groups$it
  grec <- json$groups$rec

  json$order <- list(list(mp_nodes, fda_nodes), list(it_nodes, list()), list(rec_nodes, list()))

  # order <- list(list(gmp$nodes, gfda$nodes), list(list(git$nodes)), list(list(grec$nodes)))
  # json$order <- order

  groups <- list(gmp, gfda, git, grec)
  json$groups <- groups
  json
}

#' get Category 1,2 recommendations
#' if cat2 == "yes", then it's category 2, otherwise it's category 1
#' (difference is in using given cancer type, vs. all other cancer types)
#'
#' @param MP input molecular profile
#' @param cancer_type cancer type selection
#' @param drugs_PO_FDA_biomarkers drugs data frame
#' @param cat2
#'
#' @importFrom dplyr filter
get_cat_1_2 <- function(MP,
                        cancer_type,
                        drugs_PO_FDA_biomarkers,
                        cat2 = "no")
{

  # load(paste0(getwd(), "/data/database_inputs_to_app.RData"))
  # print(cancer_type)

  ##change cat2 to lower case, in case someone typed "Yes," "YES" etc
  cat2 <- tolower(cat2)
  if(cat2 == "y" | cat2 == "ye")
  {
    cat2 <- "yes"
  }

  Type1_2_df <- NULL
  if(cat2 != "yes")
  {
    Type1_2_df <- filter(drugs_PO_FDA_biomarkers,
                                Gene.Protein %in% MP$Gene_protein,
                                Disease == cancer_type)
  } else {
    Type1_2_df <- filter(drugs_PO_FDA_biomarkers,
                                Gene.Protein %in% MP$Gene_protein,
                                Disease != cancer_type)
  }

  Type1_2_df <- check_WT_biomarkers(Type1_2_df, MP)

  Type1_2_df <- Type1_2_df[,c("Drug","Gene.Protein","Type","Alteration","Disease")]
  colnames(Type1_2_df)[2] <- "Gene or Protein"
  ##change "Disease" to "Tumor in which it is approved"
  colnames(Type1_2_df)[colnames(Type1_2_df) == "Disease"] <-
    "Tumor in which it is approved"

  if(nrow(Type1_2_df) == 0)
  {
    Type1_2_df <-
      data.frame(Note = "There are no recommended therapies in this category.")
  }
  else {
    pdrugs <- lapply(Type1_2_df$Drug, function(d) {
      temp <- d
      if(d %in% names(prop2non_prop)) {
        temp <- prop2non_prop[[d]]
      }

      temp
    })

    cat("Updating names \n")
    test <- pdrugs
    names(test) <- Type1_2_df$Drug
    print(test)

    Type1_2_df$Drug <- unlist(pdrugs)

  }

  Type1_2_df
}

#' get Category 3,4 recommendations
#' if cat4 == "yes", then it's category 4, otherwise it's category 3
#' (difference is in using given cancer type, vs. all other cancer types)
#'
#' @param MP input molecular profile
#' @param cancer_type cancer type selection
#' @param drugs_PO_FDA_biomarkers drugs data frame
#' @param cat_drugs character string - whether only targeted cancer therapies, all FDA approved therapies, or all drugs in DrugBank
#' @param list_paths list of pathways
#' @param Onc_df data frame with oncogenes and cancer type
#' @param drugs_PO_FDA_biomarkers FDA-approved precision oncology drugs with listed biomarkers
#' @param drugs_PO_FDA_targets FDA-approved precision oncology drugs with targets
#' @param drug_targets data frame of targets for specific genes/proteins - right now, have DrugBank_targets
#' @param FDA_approved_drugs list of FDA-approved drugs)
#' @param Type1 results of get_cat_1_2
#' @param Type2 results of get_cat_1_2
#' @param Type3 results of get_cat_3_4 (only use in cat4 == "yes)
#' @param cat4 defaults to "no"
#' @param subtype choose subtype of connections, if have multiple ones in list_paths
#' @return data frame
#'
#' @importFrom dplyr filter
#'
get_cat_3_4 <- function(MP,
                        cancer_type,
                        cat_drugs,
                        list_paths,
                        Onc_df,
                        drugs_PO_FDA_biomarkers,
                        drugs_PO_FDA_targets,
                        drug_targets,
                        FDA_approved_drugs,
                        Type1,
                        Type2,
                        Type3 = NULL,
                        cat4 = "no",
                        subtype = "activation")
{

  # load(paste0(getwd(), "/data/database_inputs_to_app.RData"))

  ##change cat4 to lower case, in case someone typed "Yes," "YES" etc
  cat4 <- tolower(cat4)
  if(cat4 == "y" | cat4 == "ye")
  {
    cat4 <- "yes"
  }

  cancer_path <- NULL
  if(cat4 != "yes")
  {
    ##get the pathway for this cancer
    ##get just the connections in that pathway - note that only keeping activations!
    cancer_path <- get_connections_paths(list_paths, cancer_type, subtype=subtype)
  } else {
    ##get the connections from all the other cancer pathways
    other_paths <- setdiff(unique(Onc_df$Name),
                           cancer_type)

    cancer_path <- get_connections_paths(list_paths, other_paths, subtype=subtype)
  }

  if(cat4 != "yes")
  {
    ##get the oncogenes from that pathway for all cancer types)
    Onc_df <- filter(Onc_df,
                            Name == cancer_type)
    ##make sure these oncogenes are not wild type in the MP
    ##(which means they would not in fact be oncogenic)
    ##get molecular alterations that are not wild-type in the MP
    MP_not_WT <- remove_WT(MP)
  } else
  {
    ##for this category, consider all oncogenes from all cancer types
    ##make sure these oncogenes are not wild type in the MP
    ##(which means they would not in fact be oncogenic)
    ##get molecular alterations that are not wild-type in the MP
    MP_not_WT <- remove_WT(MP)
  }

  inputs <- unique(as.character(MP_not_WT$Gene_protein))

  ##get genes among the MP that are both 1) oncogenes and 2) not wild type
  ##these will be the only inputs considered further
  ##(will look downstream of them)
  inputs_in_path <- filter(Onc_df,
                                  Oncogene %in% MP_not_WT$Gene_protein)
  inputs_in_path <- unique(as.character(inputs_in_path$Oncogene))

  ##keep just the genes that are downstream from inputs + the inputs
  downstream_genes <- get_downstream_from_inputs(cancer_path,
                                                 MP,
                                                 inputs_in_path)
  keep_genes <- unique(c(downstream_genes, inputs))

  ##just keep connections to and from these genes + inputs
  cancer_path <- filter(cancer_path,
                               (to %in% keep_genes) |
                                 (from %in% keep_genes))

  ##get the drugs which either target genes/proteins in keep_genes or for which
  ##the genes/proteins in keep_genes are biomarkers
  drugs_PO_FDA_biomarkers_keep <- filter_drugs_select_genes(drugs_PO_FDA_biomarkers,
                                                            keep_genes)
  drugs_PO_FDA_targets_keep <- filter_drugs_select_genes(drugs_PO_FDA_targets,
                                                         keep_genes)
  drugs_targets_keep <- filter_drugs_select_genes(drug_targets,
                                                  keep_genes)

  ##now get the drugs based on the chosen category
  ##by default, include everything
  drugs_biomarkers_targets <-
    rbind(drugs_PO_FDA_biomarkers_keep[,1:6],
          drugs_PO_FDA_targets_keep[,1:6],
          drugs_targets_keep)
  if(cat_drugs == "Only FDA-approved targeted therapies for cancer")
  {
    drugs_biomarkers_targets <-
      rbind(drugs_PO_FDA_biomarkers_keep[,1:6],
            drugs_PO_FDA_targets_keep[,1:6])
  }
  if(cat_drugs == "Only FDA-approved therapies")
  {
    drugs_biomarkers_targets <-
      drugs_biomarkers_targets %>% filter(Drug %in% FDA_approved_drugs)
  }

  ##remove rows with wild-type biomarkers if those biomarkers are among the (non-wild-type) inputs
  drugs_biomarkers_targets <-
    drugs_biomarkers_targets %>%
    filter(!(tolower(Alteration) %in% c("wild type","wildtype","wild","wt") &
               Gene.Protein %in% inputs))

  ##add shortest path from inputs to targets/biomarkers in drugs_biomarkers_targets
  drugs_biomarkers_targets <- add_shortest_path(drugs_biomarkers_targets,
                                                inputs, inputs_in_path,
                                                cancer_path)

  ##simplify/remove some column names for display purposes
  drugs_biomarkers_targets <- pretty_drugs_mat(drugs_biomarkers_targets)

  ##create data frame of edges by combining pathway and gene-drug info
  edges_df <- data.frame(Source = c(as.character(cancer_path$from), as.character(drugs_biomarkers_targets$Drug)),
                         Target = c(as.character(cancer_path$to), as.character(drugs_biomarkers_targets$Gene.Protein)))

  ##change "Disease" to "Tumor in which it is approved"
  colnames(drugs_biomarkers_targets)[colnames(drugs_biomarkers_targets) == "Disease"] <- "Tumor in which it is approved"

  ##take out all genes/proteins already present above
  if(cat4 != "yes")
  {
    drugs_biomarkers_targets <- remove_known_recs(drugs_biomarkers_targets,
                                                  Type1_df=Type1,
                                                  Type2_df=Type2)
  } else {
    drugs_biomarkers_targets <- remove_known_recs(drugs_biomarkers_targets,
                                                  Type1_df=Type1,
                                                  Type2_df=Type2,
                                                  Type3_df=Type3)
  }

  if(nrow(drugs_biomarkers_targets) > 0) {
    pdrugs <- lapply(drugs_biomarkers_targets$Drug, function(d) {
      temp <- d
      if(d %in% names(prop2non_prop)) {
        temp <- prop2non_prop[[d]]
      }

      temp
    })
    cat("Updating names \n")
    test <- pdrugs
    names(test) <- drugs_biomarkers_targets$Drug
    print(test)

    drugs_biomarkers_targets$Drug <- unlist(pdrugs)
  }

  colnames(drugs_biomarkers_targets)[c(2,3)] <- c("Gene or Protein","Type")

  list(drugs_mat = drugs_biomarkers_targets,
       edges = edges_df)
}
