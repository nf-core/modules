COPY
(
  SELECT
    1 AS value
)
TO
'${prefix}.parquet'
(FORMAT parquet, COMPRESSION zstd)
