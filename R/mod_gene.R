
#' gene view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS
#' @importFrom reactable reactableOutput
mod_gene_ui <- function(id) {
  ns <- NS(id)

  gene_view_grid <- create_grid(
    areas = rbind(c("left", "right")),
    cols_width = c("30%", "auto"),
    rows_height = c("auto")
  )

  grid(
    gene_view_grid,
    left = reactableOutput(ns("gene_exp_table")) %>%
      shinycssloaders::withSpinner(),
    right = plotOutput(ns("gene_exp")) %>%
      shinycssloaders::withSpinner()
  )
}

#' gene view Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom forcats fct_reorder
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @import reactable
#' @import stringr
#' @noRd
mod_gene_server <- function(id, conn, select) {
  moduleServer(id, function(input, output, session) {
    anno <- reactive({
      validate(need(select, message = "Waiting selection"))
      tbl(conn, "anno") %>%
        dplyr::filter(gene_name == !!select)
    })

    output$gene_exp_table <- renderReactable({
      dge <- anno() %>%
        select(gene_id, gene_name) %>%
        distinct() %>%
        left_join(tbl(conn, "dge"), by = "gene_id") %>%
        select(contrasts, log2FoldChange, padj) %>%
        collect()

      validate(need(any(!is.na(dge$log2FoldChange)), "Gene not tested for DE."))
      dge %>%
        mutate_at(vars(padj, log2FoldChange), ~ format(round(., digits = 2), nsmall = 2)) %>%
        mutate(contrasts = fct_reorder(contrasts, nchar(contrasts))) %>%
        arrange(desc(contrasts)) %>%
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
          columns = list(
            contrasts = colDef(
              width = 200,
            )
          ),
          defaultColDef = colDef(width = 80)
        )
    })

    output$gene_exp <- renderPlot({
      dge <- conn %>%
        tbl("dge") %>%
        collect()

      gene_l2fc <- dge %>%
        filter(gene_name == !!select) %>%
        dplyr::select(gene_name, contrasts, log2FoldChange) %>%
        mutate(y = 0, yend = 0.7)

      validate(need(nrow(gene_l2fc) > 0, "Gene not tested for DE."))

      dge %>%
        mutate(contrasts = fct_reorder(contrasts, nchar(contrasts))) %>%
        ggplot(., aes(y = contrasts, x = log2FoldChange)) +
        ylab("Density") +
        geom_density_ridges(alpha = 0.5, color = NA, bandwidth = 0.083) +
        theme_ridges() +
        geom_text(
          data = gene_l2fc,
          aes(label = paste0("↓", !!select)),
          position = position_nudge(y = 0.2),
          hjust = 0,
          colour = "red",
          # angle=45,
          size = 3.5
        ) +
        xlim(c(-2, 2))
    })
  })
}
