#' Creates PSQL table
#'
#' @import DBI
#' @import RPostgres
#'
creat_db <- function() {

  dbname <- "nmd_transcriptome"
  host <- Sys.getenv("PGHOST")
  port <- Sys.getenv("PGPORT")
  password <- Sys.getenv("PGPASSWORD")
  user <- Sys.getenv("PGUSER")
  nmd_user <- Sys.getenv("NMD_PGUSER")
  nmd_password <- Sys.getenv("NMD_PGPASSWORD")

  conn <- dbConnect(
    Postgres(),
    host = host,
    port = port,
    password = password,
    user = user
  )

  dbExecute(
    conn,
    str_glue("CREATE USER {nmd_user} WITH PASSWORD '{nmd_password}';")
  )

  dbExecute(
    conn,
    str_glue("CREATE DATABASE {dbname} WITH OWNER {nmd_user};")
  )


  dbExecute(
    conn,
    str_glue("GRANT USAGE ON SCHEMA public TO {nmd_user};")
  )


  dbExecute(
    conn,
    str_glue("GRANT SELECT ON ALL TABLES IN SCHEMA public TO {nmd_user};")
  )


  dbExecute(
    conn,
    "ALTER DEFAULT PRIVILEGES IN SCHEMA public
   GRANT SELECT ON TABLES TO {nmd_user};"
  )

  dbDisconnect(conn)
}
