- transcript name for novel should get a extra id
- DTE on database
- render_gene_card
  - DGE
  - has_support
  - n_nmd
- add gene_exp to each row of reactable so everything is aligned
- transcript_name for 7SK all wrong on anno table
- CDS source picker 

?Warning: `bindFillRole()` only works on htmltools::tag() objects (e.g., div(), p(), etc.), not objects of type 'shiny.tag.list'.


docker cache setup
RENV_PATHS_CACHE_HOST=/opt/local/renv/cache

# where the cache should be mounted in the container
RENV_PATHS_CACHE_CONTAINER=/renv/cache

# run the container with the host cache mounted in the container
docker run --rm \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    -p 14618:14618 \
    R -s -e 'renv::restore(); shiny::runApp(host = "0.0.0.0", port = 14618)'
