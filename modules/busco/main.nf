process BUSCO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::busco=5.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.2--pyhdfd78af_0':
        'quay.io/biocontainers/bsuco:5.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("$meta.id/short_summary.*.txt")                        , emit: short_summary
    tuple val(meta), path("$meta.id/run_*/full_table.tsv")        , optional:true, emit: run_full_table
    tuple val(meta), path("$meta.id/run_*/short_summary.txt")     , optional:true, emit: run_short_summary
    tuple val(meta), path("$meta.id/run_*/missing_busco_list.tsv"), optional:true, emit: run_missing_busco_list
    path "versions.yml"                                                          , emit: versions

    tuple val(meta), path("${meta.id}/run_*/full_table.tsv"),    emit: tsv
    tuple val(meta), path("${meta.id}/run_*/short_summary.txt"), emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // handling lineage
    def lineage = meta.lineage ? "--lineage_dataset ${meta.lineage}" : ""
    def autolineage_type = meta.autolineage
    if (meta.autolineage=="auto") {
        autolineage_type = "--auto-lineage"
    } else if (meta.autolineage=="prokaryote") {
        autolineage_type = "--auto-lineage-prok"
    } else if (meta.autolineage=="eukaryote") {
        autolineage_type = "--auto-lineage-euk"
    } else {
        autolineage_type = ""
    }

    """
    gzip -cdf ${fasta} > __UNCOMPRESSED_FASTA_FILE__

    export NUMEXPR_MAX_THREADS=$task.cpus

    busco \\
        --in __UNCOMPRESSED_FASTA_FILE__ \\
        --mode ${meta.mode} \\
        --out $meta.id \\
        -c $task.cpus \\
        ${lineage} \\
        ${autolineage_type} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
