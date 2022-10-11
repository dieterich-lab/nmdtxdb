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
  })


  output$gene_info <- renderUI({
    validate(need(input$gene_select, "Waiting selection"))
    req(gene_info())
  })

  output$table_transcript <- renderReactable({
    validate(need(input$gene_select, "Waiting selection"))
    gene_name <- gtf()[[1, "gene_name"]]

    conn %>%
      tbl("dtu2") %>%
      dplyr::filter(gene_name == local(gene_name)) %>%
      left_join(tbl(conn, "anno"), by = c("gene_name", "transcript_id")) %>%
      select(transcript_name, contrasts, padj, log2fold) %>%
      filter(!is.na(transcript_name)) %>%
      left_join(
        tbl(conn, "gtf") %>% select("transcript_name", "transcript_biotype"),
        by = c("transcript_name")
      ) %>%
      distinct() %>%
      collect() %>%
      mutate(contrasts = str_sub(contrasts, 5)) %>%
      reactable(
        .,
        language = reactableLang(
          filterPlaceholder = "Filter"
        ),
        filterable = TRUE,
        striped = TRUE,
        defaultSorted = c("padj"),
        showPageSizeOptions = TRUE,
        defaultPageSize = 5,
        pageSizeOptions = c(5, 10, 25, 50),
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
        defaultColDef = colDef(
          sortNALast = TRUE
        ),
        columns = list(
          padj = colDef(
            format = colFormat(digits = 2),
            filterable = FALSE
          ),
          log2fold = colDef(
            name = "log2fc",
            format = colFormat(digits = 2),
            filterable = FALSE
          ),
          contrasts = colDef(
            width = 200
          )
        )
      )
   })
}
