COPY
(
  SELECT
    *
  FROM read_csv_auto('input.csv')
)
TO
'${prefix}.parquet'
(FORMAT parquet, COMPRESSION zstd, PER_THREAD_OUTPUT)
