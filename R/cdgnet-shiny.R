possible_cancer_types <- c(
  "Colorectal cancer",
  "Pancreatic cancer",
  "Hepatocellular carcinoma",
  "Glioblastoma",
  "Thyroid cancer",
  "Acute myeloid leukemia",
  "Chronic myeloid leukemia",
  "Basal cell carcinoma",
  "Melanoma",
  "Renal cell carcinoma",
  "Bladder cancer",
  "Prostate cancer",
  "Endometrial cancer",
  "Breast cancer",
  "Small cell lung cancer",
  "Non-small cell lung cancer"
)

#' Defines shiny user interface for the cdgnet applicatio
#' @import shiny
#' @import shinyjs
#'
.cdgnetUI <- fluidPage(
  shinyjs::useShinyjs(),

  titlePanel("Therapy recommendations using biological networks"),
  HTML('<span style="color:red"><strong>Warning!</strong> The following tool is for research purposes only. It is not intended for clinical care. </span>'),
  textOutput("error"),
  tags$head(tags$style(
    "#error{color: red;
       font-style: italic;
       }"
  )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("cancer_type", "Select cancer type",
                  sort(possible_cancer_types),
                  selected = "Colorectal cancer"),
      p(tags$hr()),
      fileInput(
        "file_MP",
        "Input tsv or csv file with molecular alterations",
        accept = c(".tsv", "csv")
      ),
      helpText("If a file is not uploaded, an example profile is loaded into the app."),
      tableOutput("MP_data"),
      p(tags$hr()),
      checkboxGroupInput(
        "filterSelection",
        label = h4("Filter Recommended Therapies"),
        choices = list(
          "Same Cancer Type" = "sCancerType",
          "Same Alteration" = "sAlteration",
          "FDA Approved Drugs" = "sFDADrugs",
          "FDA Approved Targeted Cancer Drugs" = "sCancerDrugs"
        ),
        selected = c("sCancerType", "sAlteration", "sFDADrugs", "sCancerDrugs")
      ),
      p(tags$hr()),
      p(tags$label("Preprint available on ", tags$a(href="https://www.biorxiv.org/content/10.1101/605261v1", "Biorxiv"))),
      p(tags$label("Data sources"),
        tags$ul(
          tags$li(tags$a(href="https://dailymed.nlm.nih.gov/dailymed/",
                         "DailyMed")),
          tags$li(tags$a(href="https://www.genome.jp/kegg/", "KEGG")),
          tags$li(tags$a(href="https://www.drugbank.ca/", "DrugBank")),
          tags$li(tags$a(href="https://pubchem.ncbi.nlm.nih.gov/", "PubChem")),
          tags$li(tags$a(href="https://www.accessdata.fda.gov/scripts/cder/daf/",
                         "Drugs@FDA"))
        )),
      p(tags$label("Source Code on GitHub"),
        tags$ul(
          tags$li(icon("github"), tags$a(href="https://github.com/SiminaB/CDGnet", "CDGnet Source Code")),
          tags$li(icon("github"), tags$a(href="https://github.com/jkanche/nfpmShinyComponent", "CDGnet R/Sankey Visualization"))
        )
      )
    ),
    mainPanel(
      fluidRow(
        column(12,
               shinyjs::hidden(
                 wellPanel(
                   id="type1_drugs_panel",
                   textOutput("Type1_explanation"),
                   DT::dataTableOutput("Type1_drugs")
                 )
               )),
        column(12,
               shinyjs::hidden(
                 wellPanel(
                   id="type2_drugs_panel",
                   textOutput("Type2_explanation"),
                   DT::dataTableOutput("Type2_drugs")
                 )
               )),
        column(12,
               shinyjs::hidden(
                 wellPanel(
                   id="type3_drugs_panel",
                   textOutput("Type3_explanation"),
                   uiOutput("Type3_viz"),
                   DT::dataTableOutput("Type3_drugs")
                   # plotOutput("Type3_network")
                 ))
        ),
        column(12,
               shinyjs::hidden(
                wellPanel(
                  id="type4_drugs_panel",
                  textOutput("Type4_explanation"),
                  uiOutput("Type4_viz"),
                  DT::dataTableOutput("Type4_drugs")
                  # plotOutput("Type4_network")
                ))
      ),

      #download button
      downloadButton("downloadTable", "Download recommended drugs table as csv")

    )#end of Fluid row
  )#end of main panel
)
)

#' router for the application to extract parameters from the url
#' @import shiny.router
router <- make_router(
  route("/", .cdgnetUI, NA)
)

#' shiny application server logic
#' @import shiny
#' @import shiny.router
#' @import shinyjs
.cdgnetServer <- function(input, output, session) {
  print(names(list_paths_KEGG))
  router(input, output, session)

  session$userData$firstrun <- 1
  flagCategory <<- "type1" #by default, category is 1. #added by kb

  observeEvent(input$filterSelection, {
    sAlt <- FALSE
    sCType <- FALSE
    sFDA <- FALSE
    sCDrug <- FALSE

    if ("sCancerType" %in% input$filterSelection) {
      sCType <- TRUE
    }
    if ("sAlteration" %in% input$filterSelection) {
      sAlt <- TRUE
    }
    if ("sFDADrugs" %in% input$filterSelection) {
      sFDA <- TRUE
    }
    if ("sCancerDrugs" %in% input$filterSelection) {
      sCDrug <- TRUE
    }

    for(id in c("type1_drugs_panel", "type2_drugs_panel", "type3_drugs_panel", "type4_drugs_panel")) {
      shinyjs::hide(id)
    }

    if (sAlt && sCType) {
      flagCategory <<- "type1"
      shinyjs::show("type1_drugs_panel")
    }

    if (sAlt && !sCType) {
      flagCategory <<- "type2"
      shinyjs::show("type2_drugs_panel")
    }

    if (!sAlt && sCType) {
      flagCategory <<- "type3"
      shinyjs::show("type3_drugs_panel")
    }

    if (!sAlt && !sCType) {
      flagCategory <<- "type4"
      shinyjs::show("type4_drugs_panel")
    }
  })

  # parse_url <- reactive({
  #   query <- parseQueryString(session$clientData$url_search)
  #   print("in parse url")
  #   print(query)
  #   query
  # })

  MP_react <- reactive({
    input_file <- input$file_MP

    if (!is.null(input_file))
    {
      if (grepl("tsv", input_file$datapath, ignore.case = TRUE))
      {
        MP <-
          read.table(
            input_file$datapath,
            sep = "\t",
            header = TRUE,
            encoding = 'UTF-8',
            fileEncoding = 'ISO8859-1'
          )
      }
      if (grepl("csv", input_file$datapath, ignore.case = TRUE))
      {
        MP <- read.csv(input_file$datapath, header = TRUE)
      }
    } else {

      if(is.null(get_query_param()$gene)) {
        MP <- data.frame(Gene_protein = c("KRAS","PIK3CA","BRCA2"),
                         Data_type = "mutation",
                         Alteration = c("G13V","G1049S","deleterious"))
      } else {
        MP <- data.frame(Gene_protein = c(get_query_param()$gene),
                         Data_type = "mutation",
                         Alteration = c("deleterious"))
      }

    }
    ##standardize gene names from MP
    if (ncol(MP) > 1)
    {
      id2id <- gene_names_standard(as.character(MP$Gene_protein))
      MP$Gene_protein <- id2id[MP$Gene_protein]
    }

    MP
  })

  observeEvent(input$cancer_type, {
    print("inside observe event")

    if (input$cancer_type != "default") {
      from_observe_event <<- input$cancer_type
    } else {
      input$cancer_type <- "Colorectal cancer"
    }
  })

  input_reactive <- reactive({
    print(input$cancer_type)
    ctype <- input$cancer_type
    error <- ""
    print(session$userData$firstrun)
    if (session$userData$firstrun == 1) {
      print("inside first run")
      if(!is.null(get_query_param()$cancer_type)) {
        ctype <- get_query_param()$cancer_type
      }
      print(ctype)

      if(!(ctype %in% possible_cancer_types)) {
        ctype <- input$cancer_type
        error <- paste("cancer type from url does not exist, choosing default: ", ctype)
      }
      updateSelectInput(session, "cancer_type", selected=ctype)
      session$userData$firstrun <<- session$userData$firstrun + 1
    }
    print(session$userData$firstrun)
    print(ctype)
    print("end of this")
    list(type=ctype, msg=error)
  })

  output$error <- renderText(input_reactive()$msg)


  output$MP_data_explanation <-
    renderText("Alterations from molecular profiling (loaded by user):")
  output$MP_data <- renderTable(MP_react())

  ##now look at various categories of FDA-approved therapies
  ##1) Are there any FDA-approved therapies in this tumor type for these alterations
  Type1 <- reactive({
    ##cancer type
    cancer_type <- input_reactive()$type
    ##MP data
    MP <- as.data.frame(MP_react())

    Type1 <- get_cat_1_2(MP,
                         cancer_type = cancer_type,
                         drugs_PO_FDA_biomarkers)

    Type1
  })

  output$Type1_explanation <-
    renderText(
      "Category 1: FDA-approved drugs for which alterations in these genes/proteins are approved biomarkers in this tumor type:"
    )
  output$Type1_drugs <- DT::renderDataTable(DT::datatable(Type1(),
                                                          extensions = 'Buttons', options = list(
                                                            dom = 'Bfrtip',
                                                            buttons = list(
                                                              list(
                                                                extend = "csvHtml5",
                                                                title = paste0("CDGnet_Type1_drugs_", Sys.Date())
                                                              )
                                                            )
                                                          ))
                                            )

  ##2) Are there any FDA-approved therapies in other tumor types for these alterations
  Type2 <- reactive({
    ##cancer type
    cancer_type <- input_reactive()$type
    ##MP data
    MP <- as.data.frame(MP_react())

    Type2 <- get_cat_1_2(MP,
                         cancer_type = cancer_type,
                         drugs_PO_FDA_biomarkers,
                         cat2 = "yes")
    Type2
  })

  output$Type2_explanation <-
    renderText(
      "Category 2: FDA-approved drugs for which alterations in these genes/proteins are approved biomarkers in other tumor types:"
    )
  output$Type2_drugs <- DT::renderDataTable(DT::datatable(Type2(),
                                                          extensions = 'Buttons', options = list(
                                                            dom = 'Bfrtip',
                                                            buttons = list(
                                                              list(
                                                                extend = "csvHtml5",
                                                                title = paste0("CDGnet_Type2_drugs_", Sys.Date())
                                                              )
                                                            )
                                                          )))

  ##3) Are there any FDA-approved drugs for which potentially relevant genes or proteins _in the relevant pathway for this cancer_ are either biomarkers or targets? What about other (non-FDA-approved) drugs - do they target any relevant genes or proteins?
  Type3 <- reactive({
    ##cancer type
    cancer_type <- input_reactive()$type

    ##MP data
    MP <- as.data.frame(MP_react())

    ##drug types
    cat_drugs <- "Only FDA-approved targeted therapies for cancer"
    if ("sFDADrugs" %in% input$filterSelection && "sCancerDrugs" %in% input$filterSelection) {
      cat_drugs <- "Only FDA-approved targeted therapies for cancer"
    }
    else if ("sFDADrugs" %in% input$filterSelection && !"sCancerDrugs" %in% input$filterSelection) {
      cat_drugs <- "Only FDA-approved therapies"
    }
    else if (!"sFDADrugs" %in% input$filterSelection && !"sCancerDrugs" %in% input$filterSelection) {
      cat_drugs <- "All drugs in DrugBank"
    }

    Type1_df <- as.data.frame(Type1())
    Type2_df <- as.data.frame(Type2())

    Type3 <- get_cat_3_4(MP, ##input data frame
                         cancer_type = cancer_type, ##character string
                         cat_drugs = cat_drugs,
                         list_paths = list_paths_KEGG, ##list of pathways
                         Onc_df, ##data frame with oncogenes and cancer type
                         drugs_PO_FDA_biomarkers, ##FDA-approved precision oncology drugs with listed biomarkers
                         drugs_PO_FDA_targets, ##FDA-approved precision oncology drugs with targets
                         drug_targets = DrugBank_targets, ##data frame of targets for specific genes/proteins - right now, have DrugBank_targets
                         FDA_approved_drugs, ##list of FDA-approved drugs)
                         Type1 = Type1_df, ##results of get_cat_1_2
                         Type2 = Type2_df, ##results of get_cat_1_2
                         Type3 = NULL, ##results of get_cat_3_4 (only use in cat4 == "yes)
                         cat4 = "no",
                         subtype = "activation")

    Type3
  })

  output$Type3_explanation <-
    renderText(
      "Category 3: FDA-approved drugs that either have as targets these alterations or genes/proteins or as biomarkers or targets other alterations or genes/proteins in the KEGG pathway corresponding to the specific cancer type, as well as additional drugs that have the above-mentioned genes/proteins as targets. Only consider the genes/proteins downstream of the genes/proteins found altered via molecular profiling that are also specified as oncogenes (in that cancer type) for the specific cancer type in KEGG."
    )

  output$Type3_drugs <-
    DT::renderDataTable({
      table <- Type3()$drugs_mat
      if(nrow(table) > 0) {
        table$`Gene or Protein` <- paste0('<a target="_blank" href="https://oncomx.org/searchview/?gene=', table$`Gene or Protein`, '">', table$`Gene or Protein`,'</a>')
      }
      DT::datatable(table, escape = FALSE,
                    extensions = 'Buttons', options = list(
                      dom = 'Bfrtip',
                      buttons = list(
                        list(
                          extend = "csvHtml5",
                          title = paste0("CDGnet_Type3_drugs_", Sys.Date())
                        )
                      )
                    ))
    })

  output$Type3_network_caption <-
    renderText(
      "Network obtained by looking at KEGG pathway for the specific cancer type. Only the initial genes/proteins and those that are downstream of the subset of altered genes that are also assumed to be oncogenes are considered. The drugs from Categories 1, 2, and 3 are highlighted."
    )

  output$Type3_network <- renderPlot({
    inputs <-
      unique(as.character(as.data.frame(MP_react())$Gene_protein))

    Type1_df <- as.data.frame(Type1())
    Type2_df <- as.data.frame(Type2())
    Type3_df <- as.data.frame(Type3()$drugs_mat)

    drugs <- combine_drugs(Type1_df, Type2_df, Type3_df)

    edges3 <- get_edges_cats_3_4(Type3())

    graph_pathways_cats_3_4(edges3, drugs, inputs)

  })

  output$Type3_viz_caption <-
    renderText(
      "The network visualization is a very simple implementation of a sankey diagram. The focus of the visualization is the flow of evidence from the users's uploaded molecular profile/alterations leading to a recommended therapy. Each edge in the sankey diagram represents the evidence and the nodes represents either targets/therapies. The sankey diagram is accompanied by a info box on the right that shows relevant information when a node/edge is selected in the network. Selecting a drug, shows  information from pubchem (currently drug structure and publications),  whereas choosing an edge displays the evidence that was used to create the edge between the nodes (KEGG Pathways etc)."
    )

  output$Type3_viz <- renderUI({
    MP <- as.data.frame(MP_react())
    Type1_df <- as.data.frame(Type1())
    Type2_df <- as.data.frame(Type2())
    Type3_df <- as.data.frame(Type3()$drugs_mat)

    # print(Type1_df)
    # print(Type2_df)
    if (nrow(Type3()$drugs_mat) > 0) {
      res <- dataParser(MP, Type1_df, Type2_df, Type3_df)

      chart <- CdgnetViz(data = res)
      chart$render_component(shiny = TRUE)
    }
  })

  ##4) Are there any FDA-approved drugs for which potentially relevant genes or proteins _in the relevant pathway for other cancers_ are biomarkers or targets?
  Type4 <- reactive({

    ##cancer type
    cancer_type <- input_reactive()$type

    ##MP data
    MP <- as.data.frame(MP_react())

    ##drug types
    cat_drugs <- "Only FDA-approved targeted therapies for cancer"
    if ( "sFDADrugs" %in% input$filterSelection && "sCancerDrugs" %in% input$filterSelection) {
      cat_drugs <- "Only FDA-approved targeted therapies for cancer"
    }
    else if ("sFDADrugs" %in% input$filterSelection && !"sCancerDrugs" %in% input$filterSelection) {
      cat_drugs <- "Only FDA-approved therapies"
    }
    else if (!"sFDADrugs" %in% input$filterSelection && !"sCancerDrugs" %in% input$filterSelection) {
      cat_drugs <- "All drugs in DrugBank"
    }

    Type1_df <- as.data.frame(Type1())
    Type2_df <- as.data.frame(Type2())
    Type3_df <- Type3()

    Type4 <- get_cat_3_4(MP, ##input data frame
                         cancer_type = cancer_type, ##character string
                         cat_drugs = cat_drugs, ##character string - whether only targeted cancer therapies, all FDA approved therapies, or all drugs in DrugBank
                         list_paths = list_paths_KEGG, ##list of pathways
                         Onc_df, ##data frame with oncogenes and cancer type
                         drugs_PO_FDA_biomarkers, ##FDA-approved precision oncology drugs with listed biomarkers
                         drugs_PO_FDA_targets, ##FDA-approved precision oncology drugs with targets
                         drug_targets = DrugBank_targets, ##data frame of targets for specific genes/proteins - right now, have DrugBank_targets
                         FDA_approved_drugs, ##list of FDA-approved drugs)
                         Type1 = Type1_df, ##results of get_cat_1_2
                         Type2 = Type2_df, ##results of get_cat_1_2
                         Type3 = Type3_df$drugs_mat, ##results of get_cat_3_4 (only use in cat4 == "yes)
                         cat4 = "yes",
                         subtype = "activation")

    Type4
  })

  output$Type4_explanation <-
    renderText(
      "Category 4: FDA-approved drugs that either have as targets these alterations or genes/proteins or as biomarkers or targets other alterations or genes/proteins in the KEGG pathways corresponding to other cancer types, as well as additional drugs that have the above-mentioned genes/proteins as targets. Only consider the genes/proteins downstream of the genes/proteins found altered via molecular profiling that are also specified as oncogenes for various cancer types in KEGG."
    )

  output$Type4_drugs <-
    DT::renderDataTable({
      table <- Type4()$drugs_mat
      if(nrow(table) > 0) {
        table$`Gene or Protein` <- paste0('<a target="_blank" href="https://oncomx.org/searchview/?gene=', table$`Gene or Protein`, '">', table$`Gene or Protein`,'</a>')
      }

      DT::datatable(table, escape = FALSE,
        extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = list(
            list(
              extend = "csvHtml5",
              title = paste0("CDGnet_Type4_drugs_", Sys.Date())
            )
          )
        ))
    })

  output$Type4_network_caption <-
    renderText(
      "Network obtained by looking at all KEGG cancer pathways. Only the initial genes/proteins and those that are downstream of the subset of altered genes that are also assumed to be oncogenes are considered. The drugs from Categories 1, 2, 3, and 4 are highlighted."
    )

  output$Type4_network <- renderPlot({
    inputs <-
      unique(as.character(as.data.frame(MP_react())$Gene_protein))

    Type1_df <- as.data.frame(Type1())
    Type2_df <- as.data.frame(Type2())
    Type3_df <- as.data.frame(Type3()$drugs_mat)
    Type4_df <- as.data.frame(Type4()$drugs_mat)

    drugs <- combine_drugs(Type1_df, Type2_df, Type3_df, Type4_df)

    edges34 <- get_edges_cats_3_4(Type3(), Type4())

    graph_pathways_cats_3_4(edges34, drugs, inputs)

  })

  output$Type4_viz <- renderUI({
    MP <- as.data.frame(MP_react())
    Type1_df <- as.data.frame(Type1())
    Type2_df <- as.data.frame(Type2())
    Type4_df <- as.data.frame(Type4()$drugs_mat)

    if (nrow(Type4()$drugs_mat) > 0) {
      res <- dataParser(MP, Type1_df, Type2_df, Type4_df)

      chart <- CdgnetViz(data = res)
      chart$render_component(shiny = TRUE)
    }
  })

  ################################## Download handler 1 # Added by KB
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  # It should write out data to that filename.
  output$downloadTable <- downloadHandler(

    filename = function() {
      paste(flagCategory, "_drugs_table_",Sys.Date(),".csv",sep="")
    },
    content = function(file) {
      if(flagCategory == "type1") {
        #write.table(Type1()$drugs_mat, file, sep="\t", row.names = FALSE,quote=FALSE)
        write.csv(Type1(), file, row.names = FALSE, quote=FALSE, sep="\t")
      } else if (flagCategory == "type2") {
        #write.table(Type2()$drugs_mat, file, sep="\t", row.names = FALSE, quote=FALSE)
        write.csv(Type2(), file, row.names = FALSE, quote=FALSE,sep="\t")
      } else if (flagCategory == "type3") {
        #write.table(Type3()$drugs_mat, file, sep="\t", row.names = FALSE, quote=FALSE)
        write.csv(Type3()$drugs_mat, file, row.names = FALSE, quote=FALSE, sep="\t")
      } else if (flagCategory == "type4") {
        #write.table(Type4()$drugs_mat, file, sep="\t", row.names = FALSE, quote=FALSE)
        write.csv(Type4()$drugs_mat, file, row.names = FALSE, quote=FALSE, sep="\t")
      }
    } # end of content
  ) # end of download handler 1
}
