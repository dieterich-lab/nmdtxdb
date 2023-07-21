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
mod_transcript_structure_server <- function(id, conn, g_id, t_id, contrast, cds, metadata) {
  moduleServer(id, function(input, output, session) {
    anno <- tbl(conn, "anno") %>% collect()

    template <- paste(
      "<b>%{text}</b><br>",
      "<i>log2fc</i>: %{x:.2f}<br>",
      "<i>padj</i>: %{y:.2f}",
      "<extra></extra>"
    )

    data_dge <- tbl(conn, "dge") %>%
      select(-c(baseMean, lfcSE, stat, pvalue, gene_name)) %>%
      filter(
        # gene_id %in% gtf_gene_id,
        padj < 0.05, abs(log2FoldChange) > 1
      ) %>%
      filter(contrasts %in% !!contrast) %>%
      mutate(logpadj = -log10(padj + 1e-10)) %>%
      collect() %>%
      left_join(metadata %>% select(contrasts, contrast_label), by = join_by(contrasts)) %>%
      left_join(anno %>% select(gene_id, ref_gene_name), by = join_by(gene_id)) %>%
      distinct()

    data_dte <- tbl(conn, "dte") %>%
      select(-c(exonBaseMean, dispersion, stat, pvalue)) %>%
      rename(log2FoldChange = log2fold) %>%
      filter(
        # transcript_id %in% gtf_transcript_id,
        padj < 0.05, abs(log2FoldChange) > 1
      ) %>%
      filter(contrasts %in% !!contrast) %>%
      mutate(logpadj = -log10(padj + 1e-10)) %>%
      collect() %>%
      left_join(anno %>% select(gene_id, ref_gene_name), by = join_by(gene_id), multiple = "first") %>%
      left_join(metadata %>% select(contrasts, contrast_label), by = join_by(contrasts)) %>%
      distinct()

    plotly_painel <- function(x) {
      x_subset <- x %>%
        filter(target) %>%
        mutate(i = n())

      plot_ly(
        x,
        x = ~log2FoldChange,
        y = ~logpadj,
        text = ~text,
        type = "scatter",
        mode = "markers",
        hovertemplate = template,
        marker = list(
          color = "rgba(0,0,0,0.2)"
        )
      ) %>%
        add_trace(
          data = x_subset,
          marker = list(color = "red")
        ) %>%
        add_annotations(
          data = x_subset,
          showarrow = TRUE,
          arrowcolor = toRGB("red"),
          arrowhead = 3,
          font = list(size = 14, color = toRGB("red"))
        ) %>%
        add_annotations(
          text = ~ unique(x$contrast_label),
          x = 0.4,
          y = .975,
          yref = "paper",
          xref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          showarrow = FALSE,
          font = list(size = 10)
        ) %>%
        layout(
          showlegend = FALSE,
          shapes = list(
            type = "rect",
            x0 = 0,
            x1 = .8,
            xref = "paper",
            y0 = -10,
            y1 = 50,
            yanchor = 1,
            yref = "paper",
            ysizemode = "pixel",
            fillcolor = toRGB("gray80"),
            line = list(color = "transparent")
          )
        ) %>%
        config(displaylogo = FALSE, modeBarButtonsToRemove = c(
          "select2d", "drawopenpath", "drawline", "drawrect", "drawcircle",
          "eraseshape", "hoverClosestCartesian", "hoverCompareCartesian",
          "lasso2d"
        ))
    }

    output$gene_counts <- renderPlotly({
      data <- data_dge %>%
        rename(text = ref_gene_name) %>%
        mutate(target = gene_id == !!g_id)

      validate(need(nrow(data %>% filter(target)) > 0, "No DGE significant calls."))

      data %>%
        group_by(contrasts) %>%
        group_map(.f = ~ plotly_painel(.x)) %>%
        subplot(nrows = 1, shareY = TRUE) %>%
        layout(
          title = list(text = "DGE"),
          xaxis = list(title = "log2fc"),
          yaxis = list(title = "padj")
        ) %>%
        toWebGL()
    })

    output$trancript_proportions <- renderPlotly({
      data <- data_dte %>%
        mutate(
          target = transcript_id %in% !!t_id,
          text = str_c(ref_gene_name, transcript_id, sep = "-")
        )
      validate(need(nrow(data %>% filter(target)) > 0, "No DTE significant calls."))

      data %>%
        group_by(contrasts) %>%
        group_map(.f = ~ plotly_painel(.x)) %>%
        subplot(nrows = 1, shareY = TRUE) %>%
        layout(
          title = list(text = "DTE"),
          xaxis = list(title = "log2fc"),
          yaxis = list(title = "padj")
        ) %>%
        toWebGL()
    })

    output$gene_structure <- renderPlot({
      gtf <- tbl(conn, "gtf") %>%
        filter(gene_id == !!g_id) %>%
        collect()

      transcript <- gtf %>%
        dplyr::filter(type == "transcript") %>%
        filter(cds_source %in% cds) %>%
        mutate(PTC = as.character(color == "#FF0000"))

      gtf <- gtf %>%
        filter(gene_id == !!g_id) %>%
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
