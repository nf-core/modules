process DUCKDB_TABLE2PARQUET {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68261e24307fdf80b9988d6fd13cf3551735e6dc0e7e38003259f7d8efa84cdb/data'
:         'community.wave.seqera.io/library/duckdb-cli:1.5.5--9c6d18d9f687a45d' }"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*.parquet"), emit: parquet
    tuple val("${task.process}"), val('duckdb'), eval("duckdb --version | sed 's/^v//; s/ .*//'"), topic: versions, emit: versions_duckdb

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ? ", ${task.ext.args}"  : ''  // read_csv options
    def args2  = task.ext.args2  ? ", ${task.ext.args2}" : ''  // Copy options
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stem   = table.name.replaceAll(/\.(gz|zst)$/, '')
    def delim  = stem.endsWith('.tsv') ? '\\t' : stem.endsWith('.csv') ? ',' : null
    if ( ! delim ) {
        error("DUCKDB_TABLE2PARQUET: cannot determine a delimiter for '${table.name}' -- expected a .csv or .tsv suffix, optionally followed by .gz or .zst")
    }
    """
    duckdb -c "
        SET threads=${task.cpus};
        SET memory_limit='${task.memory.toGiga()}GB';
        SET temp_directory='.';
        COPY (SELECT * FROM read_csv('${table}', delim='${delim}', header=true${args})) TO '${prefix}.parquet' (FORMAT PARQUET${args2});
    "
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.parquet
    """
}
