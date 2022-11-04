adv_grid <- create_grid(
  rbind(
    c("top"),
    c("bottom_left", "bottom_right")
  ),
  c("auto", "auto", "auto"),
  c("300px", "300px")
)


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
mod_transcript_structure_ui <- function(id) {
  ns <- NS(id)

  grid(
    adv_grid,
    top = plotlyOutput(ns("gene_counts")) %>% shinycssloaders::withSpinner(),
    bottom_left = plotlyOutput(ns("trancript_proportions")) %>% shinycssloaders::withSpinner(),
    bottom_right = plotOutput(ns("gene_structure")) %>% shinycssloaders::withSpinner()
  )
}

#' phase1 Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom reactable reactableOutput getReactableState
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_split
#' @import plotly
#' @import ggtranscript
#' @import ggplot2
#' @import stringr
#' @noRd
mod_transcript_structure_server <- function(id, conn, g_id, t_name, contrast) {
  moduleServer(id, function(input, output, session) {

    output$gene_counts <- renderPlotly({
      conn %>%
        tbl("gene_counts2") %>%
        filter(gene_id == !!g_id) %>%
        collect() %>%
        plot_ly(
          type = "box",
          x = ~name,
          y = ~ log10_or_max(value),
          color = ~ factor(name),
          colors = c(I("steelblue"), I("gold"), I("forestgreen"))
        ) %>%
        config(displayModeBar = FALSE) %>%
        layout(
          title = "Gene counts",
          xaxis = list(
            title = "",
            showticklabels = FALSE
          ),
          yaxis = list(title = "log10(counts)", rangemode = "tozero"),
          legend = list(orientation = "h")
        )
    })

    output$trancript_proportions <- renderPlotly({

      tbl(conn, "tx_counts2") %>%
        filter(transcript_name %in% !!t_name) %>%
        select(-c(gene_name, transcript_id)) %>%
        collect() %>%
        plot_ly(
          type = "box",
          boxpoints = "all",
          jitter = 1,
          pointpos = 0,
          y = ~transcript_name,
          x = ~usage,
          color = ~name,
          colors = c(I("steelblue"), I("gold"), I("forestgreen")),
          orientation = "h",
          opacity = 0.8
        ) %>%
        layout(
          boxmode = "group",
          title = "Transcript proportion",
          xaxis = list(title = ""),
          showlegend = FALSE
        )
    })

    output$gene_structure <- renderPlot({
      gtf <- conn %>%
        tbl("gtf") %>%
        filter(gene_id == !!g_id) %>%
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
        # geom_range(
        #   aes(
        #     fill = transcript_biotype,
        #     height = 0.25
        #   )
        # ) +
        # geom_range(
        #   data = cds,
        #   aes(
        #     fill = transcript_biotype
        #   )
        # ) +
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
