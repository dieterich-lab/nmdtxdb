Duration: 1m 54.1s

❯ checking whether package ‘nmdtx’ can be installed ... WARNING
  See below...

❯ checking R files for non-ASCII characters ... WARNING
  Found the following file with non-ASCII characters:
    mod_gene.R
  Portable packages must use only ASCII characters in their R code,
  except perhaps in comments.
  Use \uxxxx escapes for other characters.

❯ checking Rd \usage sections ... WARNING
  Undocumented arguments in documentation object 'create_grid'
    ‘areas’ ‘cols_width’ ‘rows_height’

  Undocumented arguments in documentation object 'name'
    ‘base_path’

  Functions with \usage entries need to have the appropriate \alias
  entries, and all their arguments documented.
  The \usage entries must correspond to syntactically valid R code.
  See chapter ‘Writing R documentation files’ in the ‘Writing R
  Extensions’ manual.

❯ checking package dependencies ... NOTE
  Imports includes 24 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking for hidden files and directories ... NOTE
  Found the following hidden files and directories:
    .devcontainer
    .vscode
  These were most likely included in error. See section ‘Package
  structure’ in the ‘Writing R Extensions’ manual.

❯ checking R code for possible problems ... NOTE
  app_server: no visible binding for global variable ‘transcript_name’
  app_server: no visible binding for global variable ‘contrasts’
  app_server: no visible binding for global variable ‘padj’
  app_server: no visible binding for global variable ‘log2fold’
  app_server: no visible binding for global variable ‘.’
  mod_gene_server : <anonymous>: no visible binding for global variable
    ‘gene_name’
  mod_gene_server : <anonymous>: no visible binding for global variable
    ‘gene_id’
  mod_gene_server : <anonymous>: no visible binding for global variable
    ‘contrasts’
  mod_gene_server : <anonymous>: no visible binding for global variable
    ‘log2FoldChange’
  mod_gene_server : <anonymous>: no visible binding for global variable
    ‘padj’
  mod_gene_server : <anonymous>: no visible binding for global variable
    ‘.’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘gene_name’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘contrasts’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘transcript_name’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘padj’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘log2fold_SMG6kd_SMG7ko_control’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘log2fold_SMG5kd_SMG7ko_control’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘.’
  mod_phase1_server : <anonymous>: no visible global function definition
    for ‘colorRamp’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘transcript_id’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘value’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘total’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘type’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘start’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘end’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘transcript_biotype’
  mod_phase1_server : <anonymous>: no visible binding for global variable
    ‘strand’
  name: no visible global function definition for ‘read.csv’
  name: no visible binding for global variable ‘.’
  name: no visible binding for global variable ‘featureID’
  name: no visible binding for global variable ‘groupID’
  name: no visible binding for global variable ‘transcript_id’
  name: no visible binding for global variable ‘type’
  name: no visible binding for global variable ‘gene_id’
  name: no visible binding for global variable ‘transcript_name’
  name: no visible binding for global variable ‘gene_name’
  name: no visible binding for global variable ‘gene’
  name: no visible binding for global variable ‘SYMBOL’
  name: no visible global function definition for ‘%<>%’
  name: no visible binding for global variable ‘Replicate’
  name: no visible binding for global variable ‘Condition’
  name: no visible binding for global variable ‘Cell_line’
  name: no visible binding for global variable ‘CCG_Sample_ID’
  name: no visible binding for global variable ‘group’
  name: no visible binding for global variable ‘feature_id’
  plot_annotation: no visible binding for global variable
    ‘transcript_name’
  plot_annotation: no visible binding for global variable
    ‘transcript_biotype’
  plot_annotation: no visible binding for global variable ‘type’
  plot_annotation: no visible binding for global variable ‘start’
  plot_annotation: no visible binding for global variable ‘end’
  plot_annotation: no visible binding for global variable ‘strand’
  render_gene_card: no visible binding for global variable ‘type’
  render_gene_card: no visible binding for global variable
    ‘log2FoldChange’
  render_gene_card: no visible binding for global variable
    ‘transcript_biotype’
  render_gene_card: no visible binding for global variable ‘.’
  Undefined global functions or variables:
    %<>% . CCG_Sample_ID Cell_line Condition Replicate SYMBOL colorRamp
    contrasts end featureID feature_id gene gene_id gene_name group
    groupID log2FoldChange log2fold log2fold_SMG5kd_SMG7ko_control
    log2fold_SMG6kd_SMG7ko_control padj read.csv start strand total
    transcript_biotype transcript_id transcript_name type value
  Consider adding
    importFrom("grDevices", "colorRamp")
    importFrom("stats", "contrasts", "end", "start")
    importFrom("utils", "read.csv")
  to your NAMESPACE file.

0 errors ✔ | 3 warnings ✖ | 3 notes ✖
