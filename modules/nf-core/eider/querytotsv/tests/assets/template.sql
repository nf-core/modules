COPY
(
  SELECT
    1 AS value
)
TO
'${prefix}.tsv'
(HEADER, DELIMITER '\t')
