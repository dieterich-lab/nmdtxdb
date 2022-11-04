INITIAL_CONTRAST <- c(
  "HEK_SMG7-KO_SMG5-KD_Z245",
  "HEK_SMG7-KO_SMG6-KD_Z319"
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

  gene_info <- reactiveVal()

  updateSelectizeInput(
    session,
    "gene_select",
    choices = tbl(conn, "gtf") %>% pull("gene_name") %>% unique() %>% sort(),
    server = TRUE,
    selected = ""
  )

  updateSelectizeInput(
    session,
    "contrast_select",
    choices = tbl(conn, "metadata") %>%
      pull("group") %>%
      unique() %>%
      grep('Luc', ., value = TRUE, invert = TRUE) %>%
      sort(),
    server = TRUE,
    selected = INITIAL_CONTRAST,
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

  observeEvent(
    c(anno(), contrast()),
    ignoreNULL = TRUE,
    {

      #send_toast(msg = "Loading selection.", error = "warning", session = session)

      gene_id <- anno()[[1, "gene_id"]]
      gene_name <- anno()[[1, "gene_name"]]
      transcript_id <- anno()[["transcript_id"]]
      transcript_name <- anno()[["transcript_name"]]
      gene_info(render_gene_card(gene_id, conn))
      mod_gene_server("mod_gene1", conn, gene_name, contrast())
      mod_transcript_structure_server("mod_transcript_structure", conn, gene_name, transcript_name, contrast())
      mod_transcript_server("mod_transcript1", conn, transcript_name, contrast())
    }
  )

  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })

  # observeEvent(input$contrast_select, {
  #   contrast <- input$contrast_select
  # })


  # mod_transcript_structure_server(
  #   "mod_transcript_structure", conn, gene_name
  # )
  # mod_transcript_server(
  #   "mod_transcript1", conn, gene_name
  # )
}
