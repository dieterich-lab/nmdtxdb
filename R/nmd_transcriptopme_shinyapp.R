options(shiny.autoreload = TRUE, warn = 0)
suppressPackageStartupMessages({
  library(shiny)
  library(promises)
  library(future)
})
source("nmdtx_utils.r")

library(dplyr)
library(shiny.semantic)
library(ggplot2)
library(reactable)
library(plotly)
library(crosstalk)
library(ggtranscript)
library(stringr)


plan(multisession)

icon_with_text <- function(text = "") {
  div(icon("question circle")) %>%
    htmltools::tagAppendAttributes(., "data-tooltip" = as.character(text))
}


conn <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "nmd_transcriptome",
  host = "***REMOVED***",
  port = ***REMOVED***,
  password = "***REMOVED***",
  user = "***REMOVED***"
)

log10_or_max <- function(x) {
  log10x <- log10(x)
  log10x[log10x <= -Inf] <- min(log10x[is.finite(log10x)])
  log10x[log10x >= Inf] <- max(log10x[is.finite(log10x)])
  log10x
}


send_toast <- function(msg, session, position = "top right", class = "warning", icon = "exclamation") {
  toast(
    msg,
    class = class,
    session = session,
    toast_tags = list(position = position, showIcon = icon)
  )
}

plot_coverage <- function(gtf) {
  region <- GenomicRanges::reduce(GenomicRanges::GRanges(gtf))
  Signac::BigwigTrack(region, bigwig) +
    theme_minimal()
  # lims(x = layer_scales(p1)$x$get_limits())
}

plot_annotation <- function(gtf) {
  gtf <- gtf %>%
    mutate(name = paste0(transcript_name, " ", transcript_biotype))

  exons <- gtf %>% dplyr::filter(type == "exon")
  cds <- gtf %>% dplyr::filter(type == "CDS")
  introns <- ggtranscript::to_intron(exons, group_var = "transcript_id")

  p1 <- exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = name
    )) +
    ggtranscript::geom_range(
      aes(
        height = 0.25
      )
    ) +
    ggtranscript::geom_range(
      data = cds
    ) +
    ggtranscript::geom_intron(
      data = introns,
      aes(strand = strand)
    ) +
    theme_minimal() +
    labs(y = "") +
    theme(
      axis.ticks = element_blank(),
      legend.position = c(0.87, 0.75)
    )
}

get_gene_info <- function(gene_id) {
  url <- httr::modify_url(
    url = "https://mygene.info",
    path = paste0("v3/gene/", gene_id),
    query = "fields=symbol,summary,alias,uniprot,name"
  )

  resp <- httr::GET(url, httr::accept("application/json"))
  jsonlite::fromJSON(httr::content(resp, "text"), simplifyVector = TRUE)
}
#
# if (http_type(resp) != "application/json") {
#   stop("API did not return json", call. = FALSE)
# }
# if (http_status(resp)$category  != 'Success'){
#   stop("API request failed", call. = FALSE)
# }

adv_grid <- grid_template(
  default = list(
    areas = rbind(
      c("top_left", "top_left", "top_right"),
      c("bottom_left", "bottom_mid", "bottom_right")
    ),
    rows_height = c("50%", "50%"),
    cols_width = c("20%", "30%", "50%")
  )
)

bigwig <- list(
  control1 = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/hg38/33F-1.bigWig",
  SMG6kd_SMG7ko = "https://trackhub.dieterichlab.org/tbrittoborges/nmd_transcriptome/hg38/33F-14.bigWig"
)

render_gene_card <- function(gene_id) {
  parsed <- get_gene_info(gene_id)
  transcripts <- tbl(conn, "gtf") %>%
    filter(type == "transcript", gene_id == !!gene_id) %>%
    collect()

  dge_l2fc <- tbl(conn, "dge") %>%
    filter(gene_id == !!gene_id) %>%
    collect() %>%
    pull(log2FoldChange)
  message(dge_l2fc)
  is_nmd <- transcripts %>% with(table(transcript_biotype == "nonsense_mediated_decay"))
  n_nmd <- ifelse(length(is_nmd) == 2, is_nmd[2], 0)
  card(
    align = "left",
    div(
      a(class = "ui blue ribbon label", "Gene information"),
      class = "content",
      br(),
      br(),
      div(class = "header", parsed$symbol),
      div(parsed$name),
      hr(),
      div(
        class = "meta",
        HTML(
          stringr::str_interp("Also known as: ${paste(parsed$alias, collapse=', ')}")
        )
      ),
      p(str_glue("NMD isoforms: {n_nmd}/{sum(is_nmd)}")),
      div(
        str_glue("Up-regulated in {sum(dge_l2fc > 1)}/{length(dge_l2fc)} cohorts"),
        icon("question circle"),
      ) %>% htmltools::tagAppendAttributes(., "data-tooltip" = "Number of datasets with l2fc > 1"),
      # p("Novel transcripts : X/XX"),
      p("Single molecule evidence: X/XX"),
      p(
        "Ensembl: ",
        a(
          gene_id,
          href = paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", gene_id),
          target = "_blank"
        )
      ),
      conditionalPanel(
        !is.null(parsed$uniprot$`Swiss-Prot`),
        p(
          "Uniprot: ",
          a(
            parsed$uniprot$`Swiss-Prot`,
            href = paste0("https://www.uniprot.org/uniprotkb/", parsed$uniprot$`Swiss-Prot`),
            target = "_blank"
          )
        )
      )
    )
  )
}

ui <- semanticPage(
  sidebar_layout(
    sidebar_panel(
      width = 3,
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
      div(
        style = "position:absolute;right:1em;",
        actionButton("help", label = "", icon = icon("question"), )
      ),
      tabset(
        tabs = list(
          list(menu = "Gene expression", content = reactableOutput("gene_view") %>% shinycssloaders::withSpinner(), id = "gene_view"),
          list(menu = "Transcript table", content = plotOutput("annotation") %>% shinycssloaders::withSpinner(), id = "tx_tab"),
          list(menu = "Advanced view", content = mod_phase1_ui("mod_phase1") %>% shinycssloaders::withSpinner(), id = "advanced_tab")
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


  gtf <- reactiveVal()
  gene_info <- reactiveVal()
  coverage <- reactiveVal()
  annotation <- reactiveVal()

  gtf <- eventReactive(
    input$gene_select,
    {
      if (is.na(input$gene_select) | input$gene_select == "") {
        return(NULL)
      }
      mod_phase1_server("mod_phase1", conn, input$gene_select)
      tbl(conn, "gtf") %>%
        filter(gene_name == local(input$gene_select)) %>%
        collect()
    }
  )

  observeEvent(gtf(), ignoreNULL = TRUE, ignoreInit = TRUE, {
    send_toast(msg = "Loading selection.", class = "warning", session = session)
    gene_id <- gtf()[[1, "gene_id"]]
    gene_info(render_gene_card(gene_id))
    annotation(NULL)
    coverage(NULL)
    NULL
  })

  output$gene_view <- renderReactable({
    validate(need(input$gene_select, "\nWaiting selection"))
    tbl(conn, "anno") %>%
      select(gene_id, gene_name) %>%
      dplyr::filter(gene_id == local(gtf()[[1, "gene_id"]])) %>%
      distinct() %>%
      left_join(tbl(conn, "dge")) %>%
      select(contrasts, log2FoldChange, padj) %>%
      collect() %>%
      mutate_at(vars(padj, log2FoldChange), ~ format(round(., digits = 2), nsmall = 2)) %>%
      reactable(
        .,
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
        defaultColDef = colDef(width = 160)
      )
  })


  observeEvent(gene_info(), ignoreNULL = TRUE, ignoreInit = TRUE, {
    gtf <- gtf()
    future({
      annotation(NULL)
      plot_annotation(gtf)
    }) %...>%
      annotation()

    future({
      coverage(NULL)
      Sys.sleep(5)
      plot_annotation(gtf)
    }) %...>% coverage()

    NULL
  })

  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })

  output$coverage <- renderPlot({
    req(coverage())
  })

  output$annotation <- renderPlot({
    req(annotation())
  })


}

#' phase1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList fluidRow
#' @importFrom crosstalk filter_select filter_slider
#' @importFrom plotly plotlyOutput
#' @importFrom shinycssloaders withSpinner
mod_phase1_ui <- function(id) {
  ns <- NS(id)

  grid(
    adv_grid,
    top_left = reactableOutput(ns("table_transcript")),
    top_right = plotlyOutput(ns("dtu_volcano")),
    bottom_left = plotlyOutput(ns("gene_counts")),
    bottom_mid = plotlyOutput(ns("trancript_proportions")),
    bottom_right = plotOutput(ns("gene_structure"))
  )
}

#' phase1 Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom reactable reactableOutput getReactableState
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_split
#' @import ggtranscript
#' @import ggplot2
#' @import stringr
#' @noRd
mod_phase1_server <- function(id, conn, select) {
  moduleServer(id, function(input, output, session) {

    anno <- reactive({
      validate(need(select, message = "Waiting selection"))
      tbl(conn, "anno") %>%
        dplyr::filter(gene_name == !!select)
    })

    output$table_transcript <- renderReactable({
      anno() %>%
        select(gene_id, gene_name) %>%
        distinct() %>%
        left_join(tbl(conn, "dtu")) %>%
        dplyr::select(
          genomicData,
          contrasts,
          gene_id,
          transcript_id,
          padj,
          log2fold_SMG6kd_SMG7ko_control,
          log2fold_SMG5kd_SMG7ko_control
        ) %>%
        # mutate(is_nmd = transcript_id %in% is_nmd) %>%
        collect() %>%
        mutate(genomicData = str_sub(genomicData, 1, -3))  %>%
      reactable(
        .,
        defaultPageSize = 9,
        compact = TRUE,
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
        defaultColDef = colDef(width = 80),
        columns = list(
          contrasts = colDef(
            show = FALSE
          ),
          transcript_name = colDef(
            filterable = TRUE,
            width = 160
          ),
          transcript_id = colDef(
            show = FALSE
          ),
          gene_name = colDef(
            filterable = TRUE,
            width = 120
          ),
          log2fold_SMG6kd_SMG7ko_control = colDef(
            name = "l2fc_SMG67KD",
            format = colFormat(digits = 3)
          ),
          log2fold_SMG5kd_SMG7ko_control = colDef(
            name = "l2fc_SMG57KD",
            format = colFormat(digits = 3)
          ),
          padj = colDef(
            format = colFormat(digits = 3)
          ),
          gene_id = colDef(
            name = "ensembl",
            html = TRUE,
            width = 140,
            cell = JS("function(cellInfo) {
              const url = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
          ),
          genomicData = colDef(
            name = "ucsc",
            html = TRUE,
            width = 160,
            cell = JS("function(cellInfo) {
              const url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }")
          )
        )
      )
    })

    output$dtu_volcano <- renderPlotly({
      tbl(conn, "dtu") %>%
        collect() %>%
        plot_ly(
          .,
          # color = ~gene_id == unique(anno()$gene_id), #
          colors = c("gray", "red"),
          opacity = 0.4,
          x = ~log2fold_SMG5kd_SMG7ko_control,
          y = ~ -log10_or_max(padj),
          text = ~ paste0(
            "<b>", gene_name, "</b>",
            "<br><i>transcript_name</i>: ", transcript_name
          ),
          hoverinfo = "text"
        ) %>%
        add_markers() %>%
        layout(
          showlegend = FALSE,
          title = "DTU volcano",
          hoverlabel = list(align = "left")
        ) %>%
        toWebGL()
    })

    output$gene_counts <- renderPlotly({

      anno() %>%
        select(gene_id, gene_name) %>%
        distinct() %>%
        left_join(tbl(conn, "gene_counts")) %>%
        collect() %>%
        tidyr::pivot_longer(-c(gene_id, gene_name)) %>%
        mutate(group = str_sub(name, start = 1, end = -3)) %>%
        plot_ly(
          type = "box",
          x = ~group,
          y = ~ log10_or_max(value),
          color = ~ factor(group),
          colors = c(I("steelblue"), I("gold"), I("forestgreen"))
        ) %>%
        config(displayModeBar = FALSE) %>%
        layout(
          title = "Gene counts",
          xaxis = list(
            title = "",
            showticklabels = FALSE
          ),
          yaxis = list(title = "log10(counts)", rangemode='tozero'),
          legend = list(orientation = "h")
        )
    })


    output$trancript_proportions <- renderPlotly({
      anno() %>%
        left_join(tbl(conn, "tx_counts")) %>%
        select(-c(gene_name, transcript_id, gene_id)) %>%
        tidyr::pivot_longer(-c(transcript_name)) %>%
        collect() %>%
        mutate(group = str_sub(name, start = 1, end = -3)) %>%
        group_by(name) %>%
        mutate(total = sum(value, na.rm=TRUE)) %>%
        filter(total != 0) %>%
        ungroup() %>%
        mutate(usage = value / total) %>%
        collect() %>%
        plot_ly(
          type = "box",
          boxpoints = "all",
          jitter = 1,
          pointpos = 0,
          y = ~transcript_name,
          x = ~usage,
          color = ~group,
          colors = c(I("steelblue"), I("gold"), I("forestgreen")),
          orientation = 'h',
          opacity = 0.8)  %>%
        layout(boxmode = "group",
               title = "Transcript proportion",
               xaxis = list(title = ""),
               showlegend = FALSE)

    })

    output$gene_structure <- renderPlot({

      gtf <- anno() %>%
        left_join(tbl(conn, "gtf")) %>%
        collect()

      exons <- gtf %>% dplyr::filter(type == "exon")
      cds <- gtf %>% dplyr::filter(type == "CDS")
      introns <- to_intron(exons, group_var = "transcript_id")

      exons %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = transcript_name
        )) +
        geom_range(
          aes(
            fill = transcript_biotype,
            height = 0.25
          )
        ) +
        geom_range(
          data = cds,
          aes(
            fill = transcript_biotype
          )
        ) +
        geom_intron(
          data = introns,
          aes(strand = strand),
        ) +
        theme_minimal() +
        labs(y = "") +
        theme(
          axis.ticks = element_blank()
        )
    })
  })
}

shinyApp(ui = ui, server = server)
