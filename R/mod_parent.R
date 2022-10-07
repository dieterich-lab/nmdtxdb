#' #' parent UI Function
#' #'
#' #' @description A shiny Module.
#' #'
#' #' @param id,input,output,session Internal parameters for {shiny}.
#' #'
#' #' @noRd
#' #'
#' #' @importFrom shiny NS tagList
#' mod_parent_ui <- function(id){
#'   ns <- NS(id)
#'   tagList(
#'
#'   sidebar_layout(
#'     sidebar_panel(
#'       width = 2,
#'       selectizeInput(
#'         inputId = "gene_select",
#'         label = h2("Select a gene symbol"),
#'         choices = NULL,
#'         selected = NULL,
#'         multiple = FALSE,
#'         options = list(create = FALSE)
#'       ),
#'       uiOutput("gene_info")
#'     ),
#'     main_panel(
#'       width = 10,
#'       div(style = "position:absolute;right:1em;",
#'           actionButton("help", label = "", icon = icon("question"),
#'       )),
#'       tabset(
#'         tabs = list(
#'           list(menu = "NMD-score", content = textOutput("text") %>%  shinycssloaders::withSpinner(), id = "score_tab"),
#'           list(menu = "Transcript table", content = plotOutput("annotation") %>% shinycssloaders::withSpinner(), id = "tx_tab"),
#'           list(menu = "Advanced view", content = plotOutput("coverage") %>% shinycssloaders::withSpinner(), id = "advanced_tab")
#'         ),
#'         active = "second_tab",
#'         id = "transcript_tabset"
#'       )
#'
#'   )
#' )
#'
#'
#' options(shiny.autoreload = TRUE)
#' suppressPackageStartupMessages({
#'   library(shiny)
#'   # library(DBI)
#'   # library(RPostgres)
#'   # library(GenomicRanges)
#'   # library(ggtranscript)
#'   # library(patchwork)
#'   # library(Signac)
#'   # library(httr)
#'   library(promises)
#'   library(future)
#' })
#'
#' future({
#'   library(dplyr)
#'   library(shiny.semantic)
#'   library(ggplot2)
#' })
#'
#' plan(multisession)
#'
#' conn <- DBI::dbConnect(
#'   RPostgres::Postgres(),
#'   dbname = "nmd_transcriptome",
#'   host = "***REMOVED***",
#'   port = ***REMOVED***,
#'   password = "***REMOVED***",
#'   user = "***REMOVED***"
#' )
#'
#'
#' send_toast <- function(msg, session, position = "top right", class = "warning", icon = "exclamation") {
#'   toast(
#'     msg,
#'     class = class,
#'     session = session,
#'     toast_tags = list(position = position, showIcon = icon)
#'   )
#' }
#'
#' plot_coverage <- function(gtf) {
#'   region <- GenomicRanges::reduce(GenomicRanges::GRanges(gtf))
#'   Signac::BigwigTrack(region, bigwig) +
#'     theme_minimal()
#'   # lims(x = layer_scales(p1)$x$get_limits())
#' }
#'
#' plot_annotation <- function(gtf) {
#'   gtf <- gtf %>%
#'     mutate(name = paste0(transcript_name, " ", transcript_biotype))
#'
#'   exons <- gtf %>% dplyr::filter(type == "exon")
#'   cds <- gtf %>% dplyr::filter(type == "CDS")
#'   introns <- ggtranscript::to_intron(exons, group_var = "transcript_id")
#'
#'   p1 <- exons %>%
#'     ggplot(aes(
#'       xstart = start,
#'       xend = end,
#'       y = name
#'     )) +
#'     ggtranscript::geom_range(
#'       aes(
#'         height = 0.25
#'       )
#'     ) +
#'     ggtranscript::geom_range(
#'       data = cds
#'     ) +
#'     ggtranscript::geom_intron(
#'       data = introns,
#'       aes(strand = strand)
#'     ) +
#'     theme_minimal() +
#'     labs(y = "") +
#'     theme(
#'       # axis.text.y = element_blank(),
#'       axis.ticks = element_blank(),
#'       legend.position = c(0.87, 0.75)
#'     )
#' }
#'
#' get_gene_info <- function(gene_id) {
#'   url <- httr::modify_url(
#'     url = "https://mygene.info",
#'     path = paste0("v3/gene/", gene_id),
#'     query = "fields=symbol,summary,alias,uniprot,name"
#'   )
#'
#'   resp <- httr::GET(url, httr::accept("application/json"))
#'   jsonlite::fromJSON(httr::content(resp, "text"), simplifyVector = TRUE)
#' }
#' #
#' # if (http_type(resp) != "application/json") {
#' #   stop("API did not return json", call. = FALSE)
#' # }
#' # if (http_status(resp)$category  != 'Success'){
#' #   stop("API request failed", call. = FALSE)
#' # }
#'
#' bigwig <- list(
#'   control1 = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/hg38/33F-1.bigWig",
#'   SMG6kd_SMG7ko = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/hg38/33F-14.bigWig"
#' )
#'
#' render_gene_card <- function(gene_id) {
#'   parsed <- get_gene_info(gene_id)
#'
#'   card(
#'     align = "left",
#'     div(
#'       a(class = "ui blue ribbon label", "Gene information"),
#'       class = "content",
#'       br(),
#'       br(),
#'       div(class = "header", parsed$symbol),
#'       div(parsed$name),
#'       hr(),
#'       div(
#'         class = "meta",
#'         HTML(
#'           stringr::str_interp("Also known as: ${paste(parsed$alias, collapse=', ')}")
#'         )
#'       ),
#'       p("Predicted NMD isoforms: X/XX"),
#'       p("Novel transcripts : X/XX"),
#'       p("Single molecule evidence: X/XX"),
#'       p(
#'         "Ensembl: ",
#'         a(
#'           gene_id,
#'           href = paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", gene_id),
#'           target = "_blank"
#'         )
#'       ),
#'       conditionalPanel(
#'         !is.null(parsed$uniprot$`Swiss-Prot`),
#'         p(
#'           "Uniprot: ",
#'           a(
#'             parsed$uniprot$`Swiss-Prot`,
#'             href = paste0("https://www.uniprot.org/uniprotkb/", parsed$uniprot$`Swiss-Prot`),
#'             target = "_blank"
#'           )
#'         )
#'       )
#'     )
#'   )
#' }
#'
#'
#' server <- function(input, output, session) {
#'   output$text <- renderText({
#'     validate(need(input$gene_select, "Waiting selection"))
#'     paste("gene", input$gene_select)
#'   })
#'
#'   # output$gene_info <- renderUI({
#'   #   validate(need(input$gene_select, "Waiting selection"))
#'   #   req(rv$gene_id)
#'   #   # render_gene_card(annotation() %>% then(~.[[2]]))
#'   #   send_toast(msg = "done gene card", session = session)
#'   #   render_gene_card(rv$gene_id)
#'   # })
#'
#'   # debounce - https://shiny.rstudio.com/reference/shiny/1.0.4/debounce.html
#'   updateSelectizeInput(
#'     session,
#'     "gene_select",
#'     selected = NULL,
#'     choices = tbl(conn, "gtf") %>% pull("gene_name") %>% unique() %>% sort(),
#'     server = TRUE
#'   )
#'
#'   gtf <- reactiveVal()
#'   gene_info <- reactiveVal()
#'   coverage <- reactiveVal()
#'   annotation <- reactiveVal()
#'
#'   gtf <- eventReactive(
#'     input$gene_select,
#'     {
#'       if (is.na(input$gene_select) | input$gene_select == "") {
#'         return(NULL)
#'       }
#'       message(input$gene_select)
#'       tbl(conn, "gtf") %>%
#'         filter(gene_name == local(input$gene_select)) %>%
#'         collect()
#'     }
#'   )
#'
#'   observeEvent(gtf(), ignoreNULL = TRUE, ignoreInit = TRUE, {
#'     send_toast(msg = "Loading selection.", class = "warning", session = session)
#'     gene_id <- gtf()[1, "gene_id"]
#'     message(gene_id)
#'     gene_info(render_gene_card(gene_id))
#'     annotation(NULL)
#'     coverage(NULL)
#'     NULL
#'   })
#'
#'   observeEvent(gene_info(), ignoreNULL = TRUE, ignoreInit = TRUE, {
#'     gtf <- gtf()
#'     future({
#'       annotation(NULL)
#'       plot_annotation(gtf)
#'     }) %...>%
#'       annotation()
#'
#'     future({
#'       #plot_coverage(gtf)
#'       coverage(NULL)
#'       Sys.sleep(5)
#'       plot_annotation(gtf)
#'
#'     }) %...>% coverage()
#'
#'     NULL
#'   })
#'
#'   output$gene_info <- renderUI({
#'     req(gene_info())
#'   })
#'
#'   output$coverage <- renderPlot({
#'     req(coverage())
#'   })
#'
#'   output$annotation <- renderPlot({
#'     req(annotation())
#'   })
#' }
#'
#' shinyApp(ui = ui, server = server)
#'
