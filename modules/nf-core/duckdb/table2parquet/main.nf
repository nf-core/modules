process DUCKDB_TABLE2PARQUET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/534f8759d8fd8d4ee81840aee69254afa080a0269be641f81f3407ccf7f832bd/data'
:         'community.wave.seqera.io/library/duckdb:1.5.5--a35d96a433829842' }"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*.parquet"), emit: parquet
    tuple val("${task.process}"), val('duckdb'), eval("python3 -c 'import duckdb; print(duckdb.__version__)'"), topic: versions, emit: versions_duckdb
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/^Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stem   = table.name.replaceAll(/\.gz$/, '')
    def delim  = stem.endsWith('.tsv') ? '\\t' : stem.endsWith('.csv') ? ',' : null
    if ( ! delim ) {
        error("DUCKDB_TABLE2PARQUET: cannot determine a delimiter for '${table.name}' -- expected a .csv or .tsv suffix, optionally followed by .gz")
    }
    """
    python3 <<PYEOF
    import duckdb

    duckdb.sql(
        "COPY (SELECT * FROM read_csv('${table}', delim='${delim}', header=true)) TO '${prefix}.parquet' (FORMAT PARQUET)"
    )
    PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.parquet
    """
}
