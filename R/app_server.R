INITIAL_CONTRAST <- c(
  "HEK_SMG7KO_SMG6KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319",
  "HEK_SMG7KO_SMG5KD-KD_Z319-vs-HEK_SMG7KO_LucKD-KD_Z319"
)

#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import shiny.semantic
#' @import crosstalk
#' @import dplyr
#' @noRd
app_server <- function(input, output, session) {
  conn <- connect_db()

  metadata <- load_metadata(conn)

  gene_info <- reactiveVal()

  updateSelectizeInput(
    session,
    "gene_select",
    choices = tbl(conn, "anno") %>% pull("ref_gene_name") %>% unique() %>% sort(),
    server = TRUE,
    selected = "SRSF1"
  )

  updateSelectizeInput(
    session,
    "contrast_select",
    choices = metadata,
    options = list(
      valueField = "contrasts",
      labelField = "Knockdown",
      render = I("{
        option: function(item, escape) {
          return '<div>'
            + '<strong>' + escape(item.Knockdown) + '</strong>'
            + '<br>'
            + '<i>Cell-line</i>: ' + escape(item.cellline)
            + '<br>'
            + '<i>Knock-out</i>: ' + escape(item.Knockout)
            + '</div>';
        }
      }")
    ),
    server = TRUE,
    selected = INITIAL_CONTRAST,
  )

  cds_source_choices = data.frame(
    value = c("ensembl", "hek293gao", "openprot", "ribotish", "transdecoder"),
    label = c("Ensembl", "Gao et al., 2015", "OpenProt", "ribotish", "TransDecoder")
  )

  updateSelectizeInput(
    session,
    "cds_source_select",
    choices = cds_source_choices,
    options = list(
      valueField = "value",
      labelField = "label"),
    server = TRUE,
    selected = "ensembl",
  )

  send_toast(
    msg = "Server is ready. Choose a gene on the sidebar.", class = "success", session = session
  )

  anno <- reactive({
    validate(need(input$gene_select, "Waiting selection"))
    send_toast(msg = "Loading selection.", class = "warning", session = session)
    tbl(conn, "anno") %>%
      dplyr::filter(gene_name == !!input$gene_select) %>%
      collect()
  })

  contrast <- eventReactive(input$contrast_select, {
    input$contrast_select
  })

  cds_source <- eventReactive(input$cds_source_select, {
    input$cds_source_select
  })

  observeEvent(
    c(anno(), contrast(), cds_source()),
    ignoreNULL = TRUE,
    {
      mod_intro_server("intro_1")
      gene_id <- anno()[[1, "gene_id"]]
      ref_gene_id <- anno()[[1, "ref_gene_id"]]
      gene_name <- anno()[[1, "ref_gene_name"]]
      transcript_id <- anno()[["transcript_id"]]
      transcript_name <- anno()[["ref_transcript_name"]]
      gene_info(render_gene_card(ref_gene_id, conn))
      mod_gene_server(
        "mod_gene1", conn, gene_name, contrast()
      )
      mod_transcript_structure_server(
        "mod_transcript_structure", conn, gene_id, transcript_id, contrast(),
        cds_source()
      )
      mod_transcript_server(
        "mod_transcript1", conn, transcript_id, contrast(), cds_source()
      )
    }
  ) %>% debounce(1000)

  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })
}
