#' Plot Boxplot with Gene Highlight
#'
#' This function generates a boxplot with a highlighted gene based on the specified contrast.
#'
#' @param data A data frame containing the necessary columns: gene_name, contrasts, log2FoldChange.
#' @param gene_name The name of the gene to highlight.
#' @param xlims minimum, maximum values for x-axis
#'
#'
#' @return A ggplot object displaying the boxplot with the highlighted gene.
#'
#' @examples
# data <- data.frame(gene_name = c("Gene1", "Gene2", "Gene2"),
#                    contrasts = c("Contrast1", "Contrast2", "Contrast2"),
#                    log2FoldChange = c(1.5, -0.8, 2.2))
# fc_boxplot(data, "Gene2", "Contrast2")
#'
#' @import dplyr
#' @import ggplot2
fc_boxplot <- function(data, index, gene_l2fc, xlims) {
  data <- data[index, ]
  gene_name <- data[[1, "gene_name"]]
  contrast <- data[[1, "contrasts"]]
  data <- mutate(data, y = 1, label = gene_name)

  gene_l2fc <- gene_l2fc %>% filter(contrasts == !!contrast)
  text_color <- ifelse(data$log2FoldChange > 0, "red", "blue")
  p <- gene_l2fc %>%
    ggplot(aes(y = contrasts, x = log2FoldChange)) +
    geom_boxplot() +
    geom_text(
      data = data,
      aes(label = label, y = y),
      position = position_nudge(y = 0.08),
      hjust = 0,
      colour = text_color,
      size = 3.5
    ) +
    geom_vline(data=data, colour = text_color, aes(xintercept=log2FoldChange)) +
    lims(x = xlims) +
    labs(y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank()
    )
  htmltools::plotTag(p, alt = "plots")
}


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

  reactableOutput(ns("gene_exp_table")) %>%
    shinycssloaders::withSpinner()
}

#' gene view Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom forcats fct_reorder
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @import reactable
#' @import stringr
#' @noRd
mod_gene_server <- function(id, conn, gene_name, contrast) {
  moduleServer(id, function(input, output, session) {
    dge <- reactive({
      conn %>%
        tbl("dge") %>%
        filter(contrasts %in% !!contrast) %>%
        select(contrasts, log2FoldChange, padj, gene_name) %>%
        collect()
    })

    gene_l2fc <- dge() %>%
      select(gene_name, contrasts, log2FoldChange)

    output$gene_exp_table <- renderReactable({
      dge_nrow <- dge() %>%
        filter(gene_name == !!gene_name) %>%
        nrow()
      validate(need(dge_nrow > 0, "Gene not tested for DE."))

      xmin <- gene_l2fc %>%
        pull(log2FoldChange) %>%
        min() * 1.1
      xmax <- gene_l2fc %>%
        pull(log2FoldChange) %>%
        max() * 1.1

      dge <- dge() %>%
        filter(gene_name == !!gene_name) %>%
        mutate(padj = padj %>% scales::scientific()) %>%
        left_join(load_metadata(conn), by = "contrasts") %>%
        select(contrasts, padj, log2FoldChange, everything())

      reactable(
        dge,
        highlight = TRUE,
        wrap = FALSE,
        details = function(index) fc_boxplot(dge, index, gene_l2fc, c(xmin, xmax)),
        theme = reactableTheme(
          borderColor = "#dfe2e5",
          stripedColor = "#f6f8fa",
          highlightColor = "#FFFFBF",
          cellPadding = "8px 12px"
        ),
        columns = list(
          gene_name = colDef(
            show = FALSE
          ),
          Knockdown = colDef(
            show = FALSE
          ),
          cellline = colDef(
            show = FALSE
          ),
          Knockout = colDef(
            show = FALSE
          ),
          name = colDef(
            show = FALSE
          ),
          contrasts = colDef(
            width = 180,
            show = TRUE,
            html = TRUE,
            cell = JS("
    function(cellInfo) {
    const kd = cellInfo.row['Knockdown']
    const cl = cellInfo.row['cellline']
    const ko = cellInfo.row['Knockout'] || 'NoKO'
        return (
          '<div>' +
          '<strong>' + kd + '</strong> <br>' +
          '<small><i>Cell-line</i>: ' + cl +
          ';<i> KO</i>: ' + ko + ' </small></div>'
        )
      }")
          ),
          padj = colDef(
            filterable = FALSE,
            show = TRUE,
            align = "right"
          ),
          log2FoldChange = colDef(
            name = "log2fc",
            show = TRUE,
            format = colFormat(digits = 2),
            filterable = FALSE
          )
        ),
        defaultColDef = colDef(
          width = 100
        )
      )
    })
  })
}
