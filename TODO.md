

Replace https://mygene.info golem_utils_server.R line 98 by the gene description below:

library(dplyr)
library(biomaRt)

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

gff_compare <- read.table(
  file.path("phaseFinal/stringtie_merge/fix_comp_ref.merged_each.fix.gtf.tmap"),
  header = 1
)  %>% 
  select(ref_gene_id) %>% 
  distinct()
  
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'description'),
  filters = 'ensembl_gene_id',
  values = gff_compare$ref_gene_id,
  uniqueRows=TRUE)

annotLookup <- annotLookup %>% 
  mutate(description = stringr::str_remove(description, " \\[[^.]*$"))
saveRDS(annotLookup, 'phaseFinal/data/gene_description.RDS')
