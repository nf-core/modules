#!/bin/bash

set -x

SQL="CREATE TABLE s AS SELECT * FROM read_csv('/dev/stdin', delim = '\t', header = true, columns = { 'name': 'VARCHAR', 'seq': 'VARCHAR' }); CREATE VIEW sequences AS SELECT name, upper(seq) AS sequence, length(sequence) AS length, 'dna' AS alphabet FROM s; COPY sequences TO '$1.sequences.parquet' (FORMAT 'parquet', COMPRESSION 'zstd', OVERWRITE_OR_IGNORE 1, ROW_GROUP_SIZE 100000);"

seqkit fx2tab "$1" | duckdb -csv :memory: "$SQL"
