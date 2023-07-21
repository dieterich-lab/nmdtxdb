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
  mod_intro_server("intro_1")

  conn <- connect_db()

  metadata <- load_metadata(conn)

  cds_source_choices <- data.frame(
    cds_source = c("ensembl", "hek293gao", "openprot", "ribotish", "transdecoder"),
    cds_source2 = c("Ensembl", "Gao et al., 2015", "OpenProt", "Zhang et al., 2017", "TransDecoder")
  )

  gene_feat <- tbl(conn, "gene_feat") %>% collect()


  gene_info <- reactiveVal()

  updateSelectizeInput(
    session,
    "gene_select",
    choices = gene_feat %>% filter(any_dge) %>% pull("ref_gene_name"),
    server = TRUE,
    selected = "SRSF1"
  )

  updateSelectizeInput(
    session,
    "contrast_select",
    choices = metadata,
    options = list(
      valueField = "contrasts",
      labelField = "label",
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

  observeEvent(input$gene_select, {
    cds_choices <- gene_feat %>%
      filter(ref_gene_name == input$gene_select) %>%
      select(cds_source) %>%
      separate_rows(., cds_source, sep = ";") %>%
      left_join(cds_source_choices, by = "cds_source") %>%
      rename(label = cds_source2)

    updateSelectizeInput(
      session,
      "cds_source_select",
      choices = cds_choices,
      options = list(
        valueField = "cds_source",
        labelField = "label"
      ),
      server = TRUE,
      selected = cds_choices[[1, "cds_source"]]
    )
  })

  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })

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
      gene_id <- anno()[[1, "gene_id"]]
      ref_gene_id <- anno()[[1, "ref_gene_id"]]
      gene_name <- anno()[[1, "ref_gene_name"]]
      transcript_id <- anno()[["transcript_id"]]
      gene_info(render_gene_card(ref_gene_id, conn))

      mod_gene_server(
        "mod_gene1", conn, gene_name, contrast()
      )

      mod_transcript_server(
        "mod_transcript1", conn, transcript_id, contrast(), cds_source()
      )
    }
  )
}
