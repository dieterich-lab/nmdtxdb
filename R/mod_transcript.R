grid <- create_grid(
  rbind(
    c("top"),
    c("bottom_left")
  ),
  c("100%", "100%"),
  c("70%", "30%")
)



#' transcript view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS
#' @importFrom reactable reactableOutput
#' @importFrom patchwork plot_layout
mod_transcript_ui <- function(id) {
  ns <- NS(id)

  grid(
    grid,
    top = reactableOutput(ns("table_transcript")) %>%
      shinycssloaders::withSpinner(),
    bottom_left = plotOutput(ns("gene_structure")) %>% shinycssloaders::withSpinner(),
    bottom_right = NULL
  )
}


transcript_plot <- function(gtf) {
  if (nrow(gtf) < 1) {
    return("No CDS source to show.")
  }

  exons <- gtf %>% dplyr::filter(type == "exon")
  cds <- gtf %>% dplyr::filter(type == "CDS")
  introns <- to_intron(exons, group_var = "Name")
  feat_colors <- c("TRUE" = "firebrick", "FALSE" = "black")

  p <- exons %>%
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
    theme_minimal()
  labs(y = "") +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "top"
    )
  return(p)
}

#' transcript view Server Functions
#' @import dplyr
#' @importFrom crosstalk SharedData
#' @importFrom purrr map
#' @importFrom reactable getReactableState renderReactable
#' @import stringr
#' @noRd
mod_transcript_server <- function(id, conn, tx, contrast, cds) {
  moduleServer(id, function(input, output, session) {
    transcripts <- tbl(conn, "gtf") %>%
      filter(transcript_id %in% !!tx, type == "transcript", cds_source %in% !!cds) %>%
      select(Name, transcript_id, cds_source, color, seqnames, start, end) %>%
      collect() %>%
      mutate(
        PTC = color == "#FF0000",
        position = position_from_gtf(.),
        trackhub_url = purrr::map(
          position,
          \(x)  create_trackhub_url(position = x)
        ) %>% unlist()
      )

    not_transcrips <- tbl(conn, "gtf") %>%
      filter(transcript_id %in% !!tx, type != "transcript") %>%
      collect() %>%
      left_join(
        x = transcripts %>% select(Name, cds_source, PTC),
        y = select(., !c(cds_source)),
        by = "Name", multiple = "all"
      )

    output$table_transcript <- renderReactable({
      dte <- tbl(conn, "dte") %>%
        filter(
          transcript_id %in% !!tx,
          contrasts %in% !!contrast
        ) %>%
        select(transcript_id, contrasts, padj, log2fold) %>%
        collect() %>%
        left_join(load_metadata(conn), by = "contrasts") %>%
        select(-name)

      l2fc <- tbl(conn, "dte") %>%
        filter(contrasts %in% !!contrast) %>%
        collect()

      anno <- tbl(conn, "anno") %>%
        filter(transcript_id %in% !!tx) %>%
        select(
          ref_gene_name, transcript_id, ref_transcript_name,
          ref_transcript_id, lr_support, class_code
        ) %>%
        collect()

      cds_position <- not_transcrips %>%
        filter(type == "CDS") %>%
        group_by(transcript_id, cds_id) %>%
        summarise(seqnames = first(seqnames), start = min(start), end = max(end)) %>%
        ungroup() %>%
        mutate(cds_position = position_from_gtf(.)) %>%
        select(transcript_id, cds_id, cds_position) %>%
        mutate(
          cds_trackhub_url = purrr::map(
            cds_position,
            \(x) create_trackhub_url(position = x)
          ) %>% unlist()
        )

      df <- anno %>%
        left_join(dte, by = "transcript_id", multiple = "all") %>%
        left_join(transcripts, by = "transcript_id", multiple = "all") %>%
        left_join(cds_position, by = "transcript_id", multiple = "all") %>%
        left_join(cds_source_choices, by = "cds_source", multiple = "all") %>%
        select(transcript_id, ref_transcript_name, PTC, lr_support, everything()) %>%
        mutate(
          log2fold = round(log2fold, 2),
          padj = padj %>% scales::scientific(),
          class_code = case_when(
            class_code == "=" ~ "Complete",
            class_code == "n" ~ "Retained introns",
            class_code %in% c("c", "j", "k") ~ "Splicing variants",
            TRUE ~ "Other"
          )
        ) %>%
        select(-c(seqnames, start, Name, end, color)) %>%
        distinct(transcript_id, cds_id, contrasts, .keep_all = TRUE)

      reactable(
        df,
        highlight = TRUE,
        wrap = FALSE,
        details = function(index) {
          if (!is.na(df[[index, "contrasts"]])) {
            data <- df[index, ]
            text_color <- ifelse(data$log2fold > 0, "red", "blue")

            p <- l2fc %>%
              filter(!is.na(log2fold), contrasts == data[["contrasts"]]) %>%
              ggplot(aes(y = contrasts, x = log2fold)) +
              geom_boxplot() +
              geom_vline(data = data, colour = text_color, aes(xintercept = log2fold)) +
              labs(y = "") +
              theme_minimal() +
              theme(
                axis.text.y = element_blank()
              )

            htmltools::plotTag(p, alt = "plots", height = 150)
          }
        },
        theme = reactableTheme(
          borderColor = "#dfe2e5",
          stripedColor = "#f6f8fa",
          highlightColor = "#FFFFBF",
          cellPadding = "8px 12px"
        ),
        defaultColDef = colDef(
          sortNALast = TRUE
        ),
        columns = list(
          cds_id = colDef(
            header = with_tooltip(
              "CDS_info", "CDS ID, provenance and URL to the trackhub."
            ),
            width = 210,
            show = TRUE,
            align = "left",
            vAlign = "center",
            html = TRUE,
            cell = JS("
function (cellInfo) {
  const cds_source = cellInfo.row['cds_source2'];
  const cds_id = cellInfo.row['cds_id'];
  const cds_trackhub_url = cellInfo.row['cds_trackhub_url'];

  if (cds_id === null) {
    return 'No CDS';
  } else {
    return `<div><a href='${cds_trackhub_url}' target='_blank'>${cds_id}</a></div><div><small><i>Source</i>: ${cds_source}</small></div>`;
  }
}")
          ),
          PTC = colDef(
            header = with_tooltip(
              "PTC", "\u2713 if the transcript has a PTC else \u274c."
            ),
            show = TRUE,
            width = 60,
            align = "center",
            vAlign = "center",
            cell = function(value) {
              if_else(value == FALSE, "\u274c", "\u2713", missing = "")
            }
          ),
          lr_support = colDef(
            header = with_tooltip(
              "LRS", "\u2713 if the transcript has long read support else \u274c."
            ),
            width = 60,
            show = TRUE,
            align = "center",
            vAlign = "center",
            cell = function(value) {
              if_else(value == "FALSE", "\u274c", "\u2713", missing = "")
            }
          ),
          transcript_id = colDef(
            header = with_tooltip(
              "transcript_id", "Identifier and link to the UCSC Genome Browser Trackhub."
            ),
            width = 160,
            show = TRUE,
            vAlign = "center",
            html = TRUE,
            cell = JS("
function (cellInfo) {
  const tid = cellInfo.row['transcript_id'];
  const url = cellInfo.row['trackhub_url'];

  if (url === undefined) {
    return tid;
  } else {
    return `<a href='${url}' target='_blank'>${tid}</a>`;
  }
}")
          ),
          trackhub_url = colDef(
            show = FALSE,
          ),
          ref_transcript_name = colDef(
            header = with_tooltip(
              "ref_transcript", "Ensembl transcript name followed by type of match between reference and assembly."
            ),
            width = 160,
            vAlign = "center",
            show = TRUE,
            html = TRUE,
            cell = JS("
function(cellInfo) {
  const cc = cellInfo.row['class_code'];
  const ti = cellInfo.row['ref_transcript_id'];
  const lab1 = '<div>' + '<strong>' + cellInfo.value + '</strong> </div>'
  const lab2 = '<div><small><i>Match</i>: ' + cc + '</small></div>';

  const url = 'http://www.ensembl.org/id/' + ti;
  return `<a href='${url}' target='_blank'>${lab1}</a>${lab2}`;

}")
          ),
          contrasts = colDef(
            header = with_tooltip(
              "contrasts", "DTE comparison in the format treatment vs control."
            ),
            width = 160,
            show = TRUE,
            html = TRUE,
            vAlign = "center",
            cell = JS("
function(cellInfo) {
  const kd = cellInfo.row['Knockdown'];
  const cl = cellInfo.row['cellline'];
  const ko = cellInfo.row['Knockout'] || 'NoKO';

  if (kd === null) {
    return 'Not tested';
  } else {
    return `<div><strong>${kd}</strong><br><small><i>Cell-line</i>: ${cl};<i> KO</i>: ${ko}</small></div>`;
  }
}")
          ),
          padj = colDef(
            header = with_tooltip(
              "padj", "DTE adjusted p-value."
            ),
            filterable = FALSE,
            show = TRUE,
            width = 100,
            align = "right"
          ),
          log2fold = colDef(
            header = with_tooltip(
              "log2fc", "DTE log2FoldChange."
            ),
            show = TRUE,
            width = 80,
            format = colFormat(digits = 2),
            filterable = FALSE
          ),
          Knockdown = colDef(
            show = FALSE
          ),
          cds_source = colDef(
            show = FALSE
          ),
          cds_source2 = colDef(
            show = FALSE
          ),
          cds_id = colDef(
            show = FALSE
          ),
          cds_position = colDef(
            show = FALSE
          ),
          cds_trackhub_url = colDef(
            show = FALSE
          ),
          Knockout = colDef(
            show = FALSE
          ),
          cellline = colDef(
            show = FALSE
          ),
          contrast_label = colDef(
            show = FALSE
          ),
          label = colDef(
            show = FALSE
          ),
          class_code = colDef(
            show = FALSE
          ),
          ref_gene_name = colDef(
            show = FALSE
          ),
          ref_transcript_name = colDef(
            show = FALSE
          ),
          ref_transcript_id = colDef(
            show = FALSE
          ),
          position = colDef(
            show = FALSE
          )
        )
      )
    })
    output$gene_structure <- renderPlot({
      transcript_plot(not_transcrips)
    })
  })
}
