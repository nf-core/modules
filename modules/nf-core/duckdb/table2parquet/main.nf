process DUCKDB_TABLE2PARQUET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4259f92248127d173a8637b5f50fc37d5dca8c34ae457f7ccebb637dac955e/data'
:         'community.wave.seqera.io/library/duckdb_python:a34b8568c8c0948e' }"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*.parquet"), emit: parquet
    tuple val("${task.process}"), val('duckdb'), eval("python3 -c 'import duckdb; print(duckdb.__version__)'"), topic: versions, emit: versions_duckdb
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/^Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ? ", ${task.ext.args}"  : ''      // read_csv options
    def args2  = task.ext.args2  ? ", ${task.ext.args2}" : ''  // Copy options
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stem   = table.name.replaceAll(/\.(gz|zst)$/, '')
    def delim  = stem.endsWith('.tsv') ? '\\t' : stem.endsWith('.csv') ? ',' : null
    if ( ! delim ) {
        error("DUCKDB_TABLE2PARQUET: cannot determine a delimiter for '${table.name}' -- expected a .csv or .tsv suffix, optionally followed by .gz or .zst")
    }
    """
    python3 <<PYEOF
    import duckdb

    con = duckdb.connect(config={'threads': ${task.cpus}, 'memory_limit': ${task.memory.toGiga()}GB, 'temp_directory': '.'})
    con.sql("COPY (SELECT * FROM read_csv('${table}', delim='${delim}', header=true${args})) TO '${prefix}.parquet' (FORMAT PARQUET${args2})")
    PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.parquet
    """
}
