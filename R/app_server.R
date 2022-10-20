#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @import shiny
#' @import dplyr
#' @import shiny.semantic
#' @import ggplot2
#' @import reactable
#' @import plotly
#' @import crosstalk
#' @import ggtranscript
#' @import stringr
#' @noRd
app_server <- function(input, output, session) {
  conn <- connect_db()
  # debounce - https://shiny.rstudio.com/reference/shiny/1.0.4/debounce.html
  updateSelectizeInput(
    session,
    "gene_select",
    selected = FALSE,
    choices = tbl(conn, "gtf") %>% pull("gene_name") %>% unique() %>% sort(),
    server = TRUE
  )
  send_toast(msg = "Server is ready. Choose a gene on the sidebar.", class = "success", session = session)

  gtf <- reactiveVal()
  gene_info <- reactiveVal()

  gtf <- eventReactive(
    input$gene_select,
    {
      if (is.na(input$gene_select) | input$gene_select == "") {
        return(NULL)
      }
      tbl(conn, "gtf") %>%
        filter(gene_name == local(input$gene_select)) %>%
        collect()
    }
  )

  observeEvent(gtf(), ignoreNULL = TRUE, ignoreInit = TRUE, {
    send_toast(msg = "Loading selection.", class = "warning", session = session)
    gene_id <- gtf()[[1, "gene_id"]]
    gene_info(render_gene_card(gene_id, conn))
    mod_gene_server("mod_gene1", conn, input$gene_select)
    mod_phase1_server("mod_phase1", conn, input$gene_select)
    mod_transcript_server("mod_transcript1", conn, input$gene_select)
  })


  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })
}
