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

  div(
    reactableOutput(ns("gene_exp_table")) %>%
      shinycssloaders::withSpinner(),
    tags$button(
      id = "gene_download",
      "Download as CSV",
      onclick = "Reactable.downloadDataCSV('mod_gene1-gene_exp_table', 'table_gene.csv')"
    )
  )
}

#' gene view Server Functions
#' @import dplyr
#' @importFrom forcats fct_reorder
#' @import reactable
#' @import stringr
#' @noRd
mod_gene_server <- function(id, db, gene_name, contrast) {
  moduleServer(id, function(input, output, session) {
    dge <- reactive({
      db[["dge"]] %>%
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
        left_join(load_metadata(db) %>% select(-contrast_label), by = "contrasts") %>%
        select(contrasts, padj, log2FoldChange, everything())

      reactable(
        dge %>% select(-name),
        highlight = TRUE,
        wrap = FALSE,
        details = function(index) {
          fc_boxplot(dge, index, gene_l2fc, c(xmin, xmax)) %>%
            htmltools::plotTag(., alt = "plots", height = 150)
        },
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
          clone = colDef(
            show = FALSE
          ),
          label = colDef(
            show = FALSE
          ),
          contrasts = colDef(
            width = 200,
            show = TRUE,
            html = TRUE,
            header = with_tooltip(
              "contrasts", "DGE comparison in the format treatment vs control."
            ),
            cell = JS("
    function(cellInfo) {
    const kd = cellInfo.row['Knockdown']
    const cl = cellInfo.row['cellline']
    const ko = cellInfo.row['Knockout'] || 'NoKO'
    const clone = cellInfo.row['clone'] || 'NA'

    return `<div><strong>${kd}</strong> <br><i>Cell line</i>: ${cl}
<br><i>Clone</i>: ${clone}; <i> KO</i>: ${ko}</small></div>`;
}")
          ),
          padj = colDef(
            filterable = FALSE,
            show = TRUE,
            align = "right",
            vAlign = "center",
            header = with_tooltip(
              "padj", "DGE adjusted p-value."
            )
          ),
          log2FoldChange = colDef(
            show = TRUE,
            format = colFormat(digits = 2),
            filterable = FALSE,
            vAlign = "center",
            header = with_tooltip(
              "log2fc", "DGE log2FoldChange."
            )
          )
        ),
        defaultColDef = colDef(
          width = 100
        )
      )
    })
  })
}
