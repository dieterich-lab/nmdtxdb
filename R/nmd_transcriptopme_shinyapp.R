options(shiny.autoreload = TRUE, warn = 2, shiny.error = recover)
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
      # axis.text.y = element_blank(),
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
    cols_width = c("33%", "33%", "34%")
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
      width = 2,
      selectizeInput(
        inputId = "gene_select",
        label = h2("Select a gene symbol"),
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
          list(menu = "Advanced view", content = mod_phase1_ui("mod_phase1"), id = "advanced_tab")
        ),
        active = "second_tab",
        id = "transcript_tabset"
      )
    )
  )
)

server <- function(input, output, session) {
  output$text <- renderText({
    validate(need(input$gene_select, "Waiting selection"))
    paste("gene", input$gene_select)
  })

  # debounce - https://shiny.rstudio.com/reference/shiny/1.0.4/debounce.html
  updateSelectizeInput(
    session,
    "gene_select",
    selected = NULL,
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
      message(input$gene_select)
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

    tbl(conn, 'anno') %>%
      select(gene_id, gene_name) %>%
      # problem:
      filter(gene_id == local(input$gene_select)) %>%
      distinct() %>%
      left_join(tbl(conn, "dge")) %>%
      select(contrasts, log2FoldChange, padj) %>%
      collect() %>%
      mutate_at(vars(padj, log2FoldChange), ~format(round(., digits=2), nsmall = 2)) %>%
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
        defaultColDef = colDef(width = 140)
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
      # plot_coverage(gtf)
      coverage(NULL)
      Sys.sleep(5)
      plot_annotation(gtf)
    }) %...>% coverage()

    NULL
  })

  output$gene_info <- renderUI({
    req(gene_info())
  })

  output$coverage <- renderPlot({
    req(coverage())
  })

  output$annotation <- renderPlot({
    req(annotation())
  })

  mod_phase1_server("mod_phase1", conn, input$gene_select)
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
    # top_mid = plotlyOutput(ns("p1")),
    top_right = plotlyOutput(ns("p2")),
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
# mod_phase1_server <- function(id, conn, res_auth) {
mod_phase1_server <- function(id, conn, gene_id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    anno <- tbl(conn, "anno")

    is_nmd <- tbl(conn, "gtf") %>%
      filter(type == "transcript", transcript_biotype == "nonsense_mediated_decay") %>%
      pull(transcript_id)

    dge <- tbl(conn, "dge") %>%
      filter(padj < 0.05) %>%
      dplyr::select(contrasts, gene_id, log2FoldChange, pvalue, padj) %>%
      left_join(anno) %>%
      collect()


    dtu <- tbl(conn, "dtu") %>%
      filter(padj < 0.05) %>%
      dplyr::select(
        genomicData,
        contrasts,
        gene_id,
        transcript_id,
        padj,
        log2fold_SMG6kd_SMG7ko_control,
        log2fold_SMG5kd_SMG7ko_control
      ) %>%
      mutate(is_nmd = transcript_id %in% is_nmd) %>%
      left_join(anno) %>%
      collect() %>%
      mutate(genomicData = str_sub(genomicData, 1, -3))

    tx_shared <- dtu %>%
      dplyr::select(gene_name, transcript_name, everything()) %>%
      SharedData$new(data = ., key = ~gene_id, group = "gene")
    gene_shared <- SharedData$new(dge, key = ~gene_id, group = "gene")

    output$table_transcript <- renderReactable({
      tx_shared %>%
        reactable(
          defaultPageSize = 9,
          compact = TRUE,
          highlight = TRUE,
          selection = "single",
          onClick = "select",
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
          defaultColDef = colDef(width = 100),
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
              cell = JS("function(cellInfo) {
              const url = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            ),
            genomicData = colDef(
              name = "ucsc",
              html = TRUE,
              cell = JS("function(cellInfo) {
              const url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }")
            )
          )
        )
    })

    output$p1 <- renderPlotly({
      plot_ly(
        gene_shared,
        x = ~log2FoldChange,
        y = ~ -log10_or_max(padj),
        text = ~ paste0(
          "<b>", gene_name, "</b>",
          "<br><i>contrasts</i>: ", contrasts
        ),
        hoverinfo = "text"
      ) %>%
        add_markers() %>%
        highlight("plotly_click", "plotly_doubleclick", color = I("red")) %>%
        layout(
          title = "DGE volcano",
          hoverlabel = list(align = "left")
        ) %>%
        toWebGL()
    })


    output$p2 <- renderPlotly({
      plot_ly(
        tx_shared,
        x = ~log2fold_SMG5kd_SMG7ko_control,
        y = ~ -log10_or_max(padj),
        text = ~ paste0(
          "<b>", gene_name, "</b>",
          "<br><i>transcript_name</i>: ", transcript_name
        ),
        hoverinfo = "text"
      ) %>%
        add_markers() %>%
        highlight("plotly_click", "plotly_doubleclick", color = I("red")) %>%
        layout(
          title = "DTU volcano",
          hoverlabel = list(align = "left")
        ) %>%
        toWebGL()
    })

    output$gene_counts <- renderPlotly({
      i <- getReactableState("table_transcript", "selected")
      validate(need(!is.na(i),
        message = "Please select an entry from the table."
      ))

      gene_id <- tx_shared$data()[[i, "gene_id"]]
      tbl(conn, "gene_counts") %>%
        dplyr::filter(gene_id == !!gene_id) %>%
        tidyr::pivot_longer(-c(gene_id)) %>%
        mutate(group = str_sub(name, start = 1, end = -3)) %>%
        collect() %>%
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
          yaxis = list(title = "log10(counts)"),
          legend = list(orientation = "h")
        )
    })


    output$trancript_proportions <- renderPlotly({
      i <- getReactableState("table_transcript", "selected")
      validate(need(!is.na(i),
        message = "Please select an entry from the table."
      ))
      gene_id <- tx_shared$data()[[i, "gene_id"]]
      tbl(conn, "tx_counts") %>%
        dplyr::filter(gene_id == !!gene_id) %>%
        tidyr::pivot_longer(-c(gene_id, transcript_id)) %>%
        mutate(group = str_sub(name, start = 1, end = -3)) %>%
        group_by(name) %>%
        mutate(total = sum(value, na.rm = TRUE)) %>%
        filter(total != 0) %>%
        ungroup() %>%
        mutate(usage = value / total) %>%
        left_join(anno) %>%
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
          orientation = "h",
          opacity = 0.8
        ) %>%
        layout(
          # boxmode = "group",
          title = "Transcript proportion",
          xaxis = list(title = ""),
          showlegend = FALSE
        )
    })

    output$gene_structure <- renderPlot({
      i <- getReactableState("table_transcript", "selected")
      validate(need(!is.na(i),
        message = "Please select an entry from the table."
      ))

      gene_id <- tx_shared$data()[[i, "gene_id"]]

      # filter for transcripts in the table
      gtf <- tbl(conn, "gtf") %>%
        filter(
          gene_id %in% local(gene_id),
          type %in% c("exon", "CDS")
        ) %>%
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
