process SEQKIT_FASTATOPARQUET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    //container "docker.io/mheuermammoth/seqkit-duckdb:latest"
    container "community.wave.seqera.io/library/duckdb-cli_seqkit:46b559b641e46efd"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0' :
    //    'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sequences.parquet"), emit: parquet
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def alphabet = "${meta.alphabet}" ?: 'dna'
    def row_group_size = 100000
    def sql = "CREATE TABLE s AS SELECT * FROM read_csv('/dev/stdin', delim = '\t', header = false, columns = { 'name': 'VARCHAR', 'seq': 'VARCHAR' }); CREATE VIEW sequences AS SELECT name, upper(seq) AS sequence, length(sequence) AS length, '${alphabet}' AS alphabet FROM s; COPY sequences TO '${prefix}.sequences.parquet' (FORMAT 'parquet', COMPRESSION 'zstd', OVERWRITE_OR_IGNORE 1, ROW_GROUP_SIZE ${row_group_size});"

    """
    seqkit \\
        fx2tab \\
        $args \\
        --threads $task.cpus \\
        $fasta \\
        | duckdb -csv :memory: "$sql"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        duckdb: \$( duckdb --version | cut -f 1 -d " " )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sql = "CREATE TABLE s ( name VARCHAR, seq VARCHAR ); CREATE VIEW sequences AS SELECT name, upper(seq) AS sequence, length(sequence) AS length, 'dna' AS alphabet FROM s; COPY sequences TO '${prefix}.sequences.parquet' (FORMAT 'parquet', COMPRESSION 'zstd', OVERWRITE_OR_IGNORE 1, ROW_GROUP_SIZE 100000);"

    """
    duckdb -csv :memory: "$sql"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        duckdb: \$( duckdb --version | cut -f 1 -d " " )
    END_VERSIONS
    """
}
