options(
  shiny.autoreload = TRUE,
  warn = 0,
  shiny.port = 3838,
  shiny.host = '0.0.0.0')


suppressPackageStartupMessages({
  library(shiny)
  library(promises)
  library(future)
})

library(dplyr)
library(shiny.semantic)
library(ggplot2)
library(reactable)
library(plotly)
library(crosstalk)
library(ggtranscript)
library(stringr)

source("mod_phase1.R")
source("mod_gene.R")
#source("mod_transcript.R")
source("nmdtx_utils.r")
plan(multisession)

ui <- semanticPage(
  sidebar_layout(
    sidebar_panel(
      width = 2,
      selectizeInput(
        inputId = "gene_select",
        label = h2("Select a gene:"),
        choices = NULL,
        selected = NULL,
        multiple = FALSE,
        options = list(create = FALSE)
      ),
      uiOutput("gene_info")
    ),
    main_panel(
      width = 10,
      # div(
      #   style = "position:absolute;right:1em;",
      #   actionButton("help", label = "", icon = icon("question"))
      # ),
      tabset(
        tabs = list(
          list(
            menu = "Gene expression",
            content = mod_gene_ui("mod_gene1"),
            id = "gene_view_tab"),
          list(
            menu = "Transcript table",
            content = reactableOutput("table_transcript") %>%
              shinycssloaders::withSpinner(),
            id = "transcript_view_tab"),
          list(
            menu = "Advanced view",
            content = mod_phase1_ui("mod_phase1"),
            id = "advanced_tab")
        ),
        active = "second_tab",
        id = "transcript_tabset"
      )
    )
  )
)

server <- function(input, output, session) {
  # debounce - https://shiny.rstudio.com/reference/shiny/1.0.4/debounce.html
  updateSelectizeInput(
    session,
    "gene_select",
    selected = FALSE,
    choices = tbl(conn, "gtf") %>% pull("gene_name") %>% unique() %>% sort(),
    server = TRUE
  )
  send_toast(msg = "Server is ready. Choose a gene on the sidebar.", class = "success", session = session)

  gtf <- reactiveVal()
  gene_info <- reactiveVal()
  # coverage <- reactiveVal()
  # annotation <- reactiveVal()

  gtf <- eventReactive(
    input$gene_select,
    {
      if (is.na(input$gene_select) | input$gene_select == "") {
        return(NULL)
      }
      tbl(conn, "gtf") %>%
        filter(gene_name == local(input$gene_select)) %>%
        collect()
    }
  )

  observeEvent(gtf(), ignoreNULL = TRUE, ignoreInit = TRUE, {
    send_toast(msg = "Loading selection.", class = "warning", session = session)
    gene_id <- gtf()[[1, "gene_id"]]
    gene_info(render_gene_card(gene_id))
    mod_gene_server("mod_gene1", conn, input$gene_select)
    mod_phase1_server("mod_phase1", conn, input$gene_select)
    #mod_transcript_server("mod_transcript1", conn, input$gene_select)
    # annotation(NULL)
    # coverage(NULL)
    NULL
  })

  # observeEvent(gene_info(), ignoreNULL = TRUE, ignoreInit = TRUE, {
  #   gtf <- gtf()
  #   future({
  #     annotation(NULL)
  #     plot_annotation(gtf)
  #   }) %...>%
  #     annotation()
  #
  #   future({
  #     coverage(NULL)
  #     Sys.sleep(5)
  #     plot_annotation(gtf)
  #   }) %...>% coverage()
  #
  #   NULL
  # })

  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })

  output$table_transcript <- renderReactable({
    validate(need(input$gene_select, "Waiting selection"))
    gene_name <- gtf()[[1, "gene_name"]]

    conn %>%
      tbl('dtu2') %>%
      dplyr::filter(gene_name == local(gene_name)) %>%
      left_join(tbl(conn, 'anno'), by=c("gene_name", "transcript_id")) %>%
      select(transcript_name, contrasts, padj, log2fold) %>%
      filter(!is.na(transcript_name)) %>%
      left_join(
        tbl(conn, "gtf") %>% select("transcript_name", "transcript_biotype"),
        by = c("transcript_name")) %>%
      distinct() %>%
      collect() %>%
      #mutate(genomicData = str_sub(genomicData, 1, -3))  %>%
      mutate(contrasts = str_sub(contrasts, 5))  %>%
      reactable(
        .,
        language = reactableLang(
          filterPlaceholder = 'Filter'),
        filterable = TRUE,
        striped = TRUE,
        defaultSorted = c("padj"),
        showPageSizeOptions = TRUE,
        defaultPageSize = 5,
        pageSizeOptions = c(5, 10, 25, 50),
        highlight = TRUE,
        wrap = FALSE,
        rowStyle = list(cursor = "pointer"),
        theme = reactableTheme(
          stripedColor = "#f6f8fa",
          highlightColor = "#f0f5f9",
          cellPadding = "8px 12px",
          rowSelectedStyle = list(
            backgroundColor = "#eee",
            boxShadow = "inset 2px 0 0 0 #FF0000"
          )
        ),
        defaultColDef = colDef(
          sortNALast = TRUE),
        columns = list(
          padj = colDef(
            format = colFormat(digits = 2)
          ),
          log2fold = colDef(
            name = 'log2fc',
            format = colFormat(digits = 2)
          ),
          contrasts = colDef(
            width = 200
          )
        #   gene_id = colDef(
        #     name = "ensembl",
        #     html = TRUE,
        #     width = 140,
        #     cell = JS("function(cellInfo) {
        #       const url = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' + cellInfo.value
        #       return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
        #     }"),
        #   ),
          # genomicData = colDef(
          #   name = "UCSC Browser",
          #   html = TRUE,
          #   width = 160,
          #   cell = JS("function(cellInfo) {
          #     const url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr' + cellInfo.value
          #     return '<a href=\"' + url + '\" target=\"_blank\">' + UCSC Browser + '</a>'
          #   }")
          # )
        )
      )
   })


  # output$gene_view <- renderPlot({
  #   req(coverage())
  # })
  #
  # output$annotation <- renderPlot({
  #   req(annotation())
  # })
}

shinyApp(ui = ui, server = server)
