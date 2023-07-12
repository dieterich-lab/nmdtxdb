adv_grid <- create_grid(
  rbind(
    c("top"),
    c("bottom_left", "bottom_right")
  ),
  c("50%", "50%", "50%"),
  c("auto", "auto", "auto")
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
    top = plotlyOutput(ns("trancript_proportions")) %>% shinycssloaders::withSpinner(),
    bottom_left = plotlyOutput(ns("gene_counts")) %>% shinycssloaders::withSpinner(),
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
mod_transcript_structure_server <- function(id, conn, g_id, t_id, contrast, cds) {
  moduleServer(id, function(input, output, session) {
    gtf <- conn %>%
      tbl("gtf") %>%
      filter(gene_id == !!g_id) %>%
      collect()

    t_id <- intersect(gtf$transcript_id, t_id)

    tx_counts <- tbl(conn, "tx_counts") %>%
      filter(transcript_id %in% !!t_id) %>%
      filter(contrasts %in% !!contrast) %>%
      collect()


    output$gene_counts <- renderPlotly({
      conn %>%
        tbl("dge") %>%
        filter(gene_id == !!g_id) %>%
        filter(contrasts %in% !!contrast) %>%
        collect() %>%
        plot_ly(
          type = "box",
          x = ~group,
          y = ~ log10_or_max(value),
          color = ~ factor(group)
        ) %>%
        config(displayModeBar = FALSE) %>%
        layout(
          title = "Gene expression",
          hovermode = FALSE,
          xaxis = list(
            title = "",
            showticklabels = FALSE
          ),
          yaxis = list(title = "log10(counts)", rangemode = "tozero"),
          legend = list(orientation = "h")
        )
    })

    output$trancript_proportions <- renderPlotly({
      plot_ly(
        tx_counts,
        type = "box",
        boxpoints = "all",
        jitter = 1,
        pointpos = 0,
        y = ~transcript_id,
        x = ~usage,
        color = ~group,
        orientation = "h"
      ) %>%
        config(displayModeBar = FALSE) %>%
        layout(
          boxmode = "group",
          hovermode = TRUE,
          # title = "Transcript usage per group",
          xaxis = list(title = ""),
          yaxis = list(title = ""),
          legend = list(orientation = "h")
        )
    })

    output$gene_structure <- renderPlot({
      transcript <- gtf %>%
        dplyr::filter(type == "transcript") %>%
        filter(cds_source %in% cds) %>%
        mutate(PTC = as.character(color == "#FF0000"))

      gtf <- gtf %>%
        select(-c(cds_source, color)) %>%
        dplyr::filter(type != "transcript") %>%
        left_join(transcript %>% select(cds_source, PTC, Name), by = "Name") %>%
        filter(!is.na(cds_source))
      validate(need(nrow(gtf) > 0, "No CDS source to show."))

      exons <- gtf %>% dplyr::filter(type == "exon")
      cds <- gtf %>% dplyr::filter(type == "CDS")
      introns <- to_intron(exons, group_var = "Name")
      feat_colors <- c("TRUE" = "firebrick", "FALSE" = "black")
      exons %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = Name
        )) +
        geom_range(
          aes(
            fill = PTC,
            height = 0.25
          )
        ) +
        geom_range(
          data = cds,
          aes(
            fill = PTC
          )
        ) +
        geom_intron(
          data = introns,
          aes(strand = strand),
        ) +
        scale_fill_manual(values = feat_colors) +
        theme_minimal() +
        facet_wrap(~cds_source, ncol = 1, scales = "free_y", drop = TRUE) +
        labs(y = "") +
        theme(
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "top"
        )
    })
  })
}
