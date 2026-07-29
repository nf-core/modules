process DUCKDB_TABLE2PARQUET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/duckdb:1.5.5--b362367d86dce804'
        : 'community.wave.seqera.io/library/duckdb:1.5.5--09cfd9fc55b2e6d4'}"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*.parquet"), emit: parquet
    tuple val("${task.process}"), val('duckdb'), eval("python3 -c 'import duckdb; print(duckdb.__version__)'"), topic: versions, emit: versions_duckdb

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
