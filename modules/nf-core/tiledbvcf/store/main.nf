
process TILEDBVCF_STORE {
    tag "$meta.id"
    label 'process_medium'

    conda "tiledb::tiledbvcf-py=0.21.0"
    container "${ 'tiledb/tiledbvcf-cli:0.21.0' }"
    containerOptions '--entrypoint ""'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta), val(tiledb_array_uri)

    output:
    tuple val(meta), path("${updated_db}")    , optional:true, emit: updateddb
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    tiledbvcf \\
        store \\
        --uri $tiledb_array_uri \\
        $vcf \\
        --threads $task.cpus \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiledbvcf: \$(tiledbvcf --version 2>&1 | head -n 1 | sed 's/^TileDB-VCF version //')
    END_VERSIONS
    """
}
