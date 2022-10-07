#' phase1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS plotOutput
#' @importFrom reactable reactableOutput
#' @importFrom plotly plotlyOutput
mod_phase1_ui <- function(id) {
  ns <- NS(id)

  grid(
    adv_grid,
    top_left = reactableOutput(ns("table_transcript")) %>% shinycssloaders::withSpinner(),
    top_right = plotlyOutput(ns("dtu_volcano")) %>% shinycssloaders::withSpinner(),
    bottom_left = plotlyOutput(ns("gene_counts")) %>% shinycssloaders::withSpinner(),
    bottom_mid = plotlyOutput(ns("trancript_proportions")) %>% shinycssloaders::withSpinner(),
    bottom_right = plotOutput(ns("gene_structure")) %>% shinycssloaders::withSpinner()
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
        left_join(tbl(conn, "dtu"), by="gene_id")  %>%
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
      gene_id <- anno() %>%
        pull(gene_id) %>%
        unique()

      tbl(conn, "dtu") %>%
        collect() %>%
        mutate(selected = ifelse(gene_name == !!select, "T", "F")) %>%
        plot_ly(
          .,
          x = ~log2fold_SMG5kd_SMG7ko_control,
          y = ~ -log10_or_max(padj),
          text = ~ paste0(
            "<b>", gene_name, "</b>",
            "<br><i>transcript_name</i>: ", transcript_name
          ),
          hoverinfo = "text"
        ) %>%
        add_markers(
          group = ~factor(selected),
          color = ~factor(selected),
          colors=colorRamp(c("gray", "red")),
          opacity = c(0.4, 0.9)) %>%
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
        left_join(tbl(conn, "gene_counts"), by='gene_id') %>%
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
        left_join(tbl(conn, "tx_counts"), by = c("transcript_id", "gene_id")) %>%
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
        left_join(tbl(conn, "gtf"), by = c("transcript_id", "gene_id", "transcript_name", "gene_name")) %>%
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
