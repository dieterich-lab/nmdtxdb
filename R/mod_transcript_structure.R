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

    metadata <- load_metadata(conn)
    metadata <- metadata %>%
      mutate(labels = str_glue_data(.,
        '<b> {Knockdown} </b> <br> <i>Cell-line</i>: {cellline}; <i>KO</i>: {Knockout}'
      ))
    meta_labeller <- setNames(as.character(metadata$labels), metadata$contrasts)

    output$gene_counts <- renderPlotly({

      data <- conn %>%
        tbl("dge") %>%
        filter(contrasts %in% !!contrast) %>%
        mutate(
          logpadj = -log10(padj + 1e-10)
        ) %>%
        collect()

      p <- data %>%
        ggplot(aes(x = log2FoldChange, y = logpadj, label = gene_name)) +
        geom_point(color = 'gray', alpha= .20) +
        geom_point(data=data %>% filter(gene_id == !!g_id), color='red') +
        theme_minimal() +
        facet_grid(
          . ~ contrasts,
          labeller =  labeller(contrasts=meta_labeller))

      p %>%
        ggplotly() %>%
        config(displaylogo = FALSE) %>%
        layout(dragmode = "select", hovermode = "x") %>%
        toWebGL()

    })

    output$trancript_proportions <- renderPlotly({
      data <- conn %>%
        tbl("dte") %>%
        filter(contrasts %in% !!contrast) %>%
        mutate(
          logpadj = -log10(padj + 1e-10)
        ) %>%
        collect()

      p <- data %>%
        ggplot(aes(x = log2fold, y = logpadj, label = gene_id)) +
        geom_point(color = 'gray', alpha= .20) +
        geom_point(data=data %>% filter(gene_id == !!g_id), color='red') +
        facet_grid(
          . ~ contrasts,
          labeller =  labeller(contrasts=meta_labeller)) +
        theme_minimal()

      p %>%
        ggplotly() %>%
        config(displaylogo = FALSE) %>%
        toWebGL()

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
