library(dplyr)
library(extrafont)

transcripts <- tbl(conn, "gtf") %>%
  filter(type == "transcript", cds_source != 'transdecoder') %>%
  mutate(PTC = color == "#FF0000")  %>%
  select(transcript_id, PTC) %>%
  group_by(transcript_id) %>%
  summarise(PTC = any(PTC)) %>%
  collect()

df <- tbl(conn, "anno") %>%
  mutate(
    class_code = case_when(
      class_code == "=" ~ "Complete",
      class_code == "n" ~ "Retained introns",
      class_code %in% c("c", "j", "k") ~ "Splicing variants",
      TRUE ~ "Other"
    )) %>%
  collect()

x <- left_join(df, transcripts, by="transcript_id") %>%
  count(class_code, PTC)

theme_Publication <- function(base_size=14, base_family = "Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}
%>%
  overview <- x %>%
  mutate(class_code = forcats::fct_reorder(class_code, desc(n))) %>%
  ggplot(aes(fill=PTC, y=n, x=class_code, label=n)) +
  geom_bar(position = "dodge", stat="identity") +
  geom_text(size = 5, position=position_dodge(0.9), color="black", vjust = -0.2) +
  scale_fill_manual(
    values = c('black', 'firebrick', 'gray50'),
    name = 'is_PTC') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                n.breaks = 5L,
                expand = c(0, 0, 0.1, 0)
  ) +
  theme_linedraw() +
  theme_Publication() +
  labs(
    title = "NMDtxDB",
    subtitle = "Dataset overview (excluding TransDecoder)",
    x = "Reference match",
    y = "Numb. of transcripts"
  )

ggsave(
  filename = "~/Nextcloud2/nmd_transcriptome_figures/overview.pdf",
  plot = overview,
  dpi = 'retina')

transcripts <- tbl(conn, "gtf") %>%
  filter(type == "transcript") %>%
  mutate(PTC = color == "#FF0000")  %>%
  select(transcript_id, PTC, cds_source) %>%
  collect()

df <- tbl(conn, "anno") %>%
  mutate(
    class_code = case_when(
      class_code == "=" ~ "Complete",
      class_code == "n" ~ "Retained introns",
      class_code %in% c("c", "j", "k") ~ "Splicing variants",
      TRUE ~ "Other"
    )) %>%
  collect()

x <- left_join(df, transcripts, by="transcript_id", multiple = "all") %>%
  count(class_code, PTC, cds_source) %>%
  distinct()

overview2 <- x %>%
  mutate(class_code = forcats::fct_reorder(class_code, desc(n))) %>%
  ggplot(aes(fill=PTC, y=n, x=class_code, label=n)) +
  geom_bar(position = "dodge", stat="identity") +
  geom_text(size = 5, position=position_dodge(0.9), color="black", vjust = -0.2) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_fill_manual(
    values = c('black', 'firebrick', 'gray50'),
    name = 'is_PTC') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                n.breaks = 5L,
                expand = c(0, 0, 0.1, 0)
  ) +
  facet_wrap(~cds_source) +
  theme_linedraw() +
  theme_Publication() +
  labs(
    title = "NMDtxDB",
    subtitle = "Dataset overview by CDS soruce",
    x = "Reference match",
    y = "Numb. of transcripts"
  )

ggsave(
  filename = "~/Nextcloud2/nmd_transcriptome_figures/overview_bysource.pdf",
  plot = overview2,
  dpi = 'retina')

ggplot(head(x), aes(fill = PTC, values = n)) +
  geom_waffle(size = .25, n_rows = 10, flip = TRUE) +
  facet_wrap(~class_code, nrow = 1, strip.position = "bottom") +
  scale_x_discrete() +
  scale_y_continuous(
    expand = c(0,0)) +
  ggthemes::scale_fill_tableau(name=NULL) +
  coord_equal()  +
  theme_minimal() + #
  theme(panel.grid = element_blank(), axis.ticks.y = element_line()) +
  guides(fill = guide_legend(reverse = TRUE))

tbl(conn, 'dge') %>%
  select(log2FoldChange, padj, contrasts) %>%
  mutate(padj = ifelse(padj <= 1e-20, 1e-20, padj)) %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point() +
  theme_linedraw() +
  theme_Publication() +
  geom_vline(xintercept=c(-0.1, 0.1), col="black", linetype='dotted') +
  geom_hline(yintercept=-log10(1e-20), col="black", linetype='dotted') +
  annotate("text", x = 0, y = -log10(1e-20) + 0.5, label = "Max -log10(prob_no_change)") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype='dotted') +
  labs(
    title = "DGE volcano",
  ) +
  facet_wrap(~contrasts, ncol = 3)




tx_anno <- tbl(conn, 'gtf') %>%
  filter(type == 'transcript', cds_source != 'transdecoder') %>%
  select(transcript_id, color) %>%
  left_join(tbl(conn, 'anno'), by="transcript_id") %>%
  group_by(gene_id, )
collect()

x <- tbl(conn, 'dge') %>%
  filter(contrasts == "HEK_NoKO_SMG5KD-KD_Z023-vs-HEK_NoKO_LucKD-KD_Z023") %>%
  collect() %>%
  left_join(tx_anno, multiple = "all", by = join_by(gene_id, gene_name)) %>%
  distinct() %>%
  mutate(PTC = color == "#FF0000") %>%
  group_by(contrasts, gene_id) %>%
  summarise(log2FoldChange = first(log2FoldChange), PTC = any(PTC), padj = first(padj))  %>%
  collect()

library(patchwork)

p2 <- ggplot(x, aes(x=log2FoldChange, fill = PTC)) +
  geom_density(alpha=0.8) +
  # facet_wrap(~contrasts, ncol = 1) +
  scale_fill_manual(
    values = c('black', 'firebrick', 'gray50'),
    name = 'is_PTC') +
  xlim(-5, 5) +
  theme_linedraw() +
  theme_Publication() +
  theme(
    plot.margin = margin(0.3, 5.5, 5.5, 5.5))

p1 <- x %>%
  mutate(padj = ifelse(padj <= 1e-20, 1e-20, padj)) %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj), fill = PTC)) +
  geom_point() +
  xlim(-5, 5) +
  scale_fill_manual(
    values = c('black', 'firebrick', 'gray50'),
    name = 'is_PTC') +
  theme_linedraw() +
  theme_Publication() +
  geom_vline(xintercept=c(-0.1, 0.1), col="black", linetype='dotted') +
  geom_hline(yintercept=-log10(1e-20), col="black", linetype='dotted') +
  annotate("text", x = 0, y = -log10(1e-20) + 0.5, label = "Max -log10(prob_no_change)") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype='dotted') +
  labs(
    title = "DGE volcano",
  ) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = margin(5.5, 5.5, -1, 5.5))

p1 / p2 + plot_layout(heights = c(1, .2))
