# observeEvent(gene_info(), ignoreNULL = TRUE, ignoreInit = TRUE, {
#   gtf <- gtf()
#   future({
#     annotation(NULL)
#     plot_annotation(gtf)
#   }) %...>%
#     annotation()
#
#   future({
#     coverage(NULL)
#     Sys.sleep(5)
#     plot_annotation(gtf)
#   }) %...>% coverage()
#
#   NULL
# })

# observeEvent(gtf(), ignoreNULL = TRUE, ignoreInit = TRUE, {
#   send_toast(msg = "Loading selection.", class = "warning", session = session)
#   gene_id <- gtf()[[1, "gene_id"]]
#   gene_info(render_gene_card(gene_id))
#   mod_gene_server("mod_gene1", conn, input$gene_select)
#   mod_transcript_structure_server("mod_transcript_structure", conn, input$gene_select)
# mod_transcript_server("mod_transcript1", conn, input$gene_select)
# annotation(NULL)
# coverage(NULL)
#   NULL
# })
