conn <- dbConnect(
  RPostgres::Postgres(),
  host = "***REMOVED***",
  port = ***REMOVED***,
  password = "WZupLe2wpejLBz7R2qrKWQRUdNvcWaGoBvcpracviYSPJMQ4duCXaGpMWHN8PjBU",
  user = "postgres"
)

dbExecute(
  conn,
  "CREATE USER ***REMOVED*** WITH PASSWORD '***REMOVED***';"
)

dbExecute(
  conn,
  "CREATE DATABASE nmd_transcriptome WITH OWNER ***REMOVED***;"
)


dbExecute(
  conn,
  "GRANT USAGE ON SCHEMA public TO ***REMOVED***;"
)


dbExecute(
  conn,
  "GRANT SELECT ON ALL TABLES IN SCHEMA public TO ***REMOVED***;"
)


dbExecute(
  conn,
  "ALTER DEFAULT PRIVILEGES IN SCHEMA public
   GRANT SELECT ON TABLES TO ***REMOVED***;"
)


dbDisconnect(conn)
