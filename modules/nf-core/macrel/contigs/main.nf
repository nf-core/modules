process MACREL_CONTIGS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::macrel=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macrel:1.1.0--py36hc5360cc_0':
        'quay.io/biocontainers/macrel:1.1.0--py36hc5360cc_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*/*.smorfs.faa.gz")      , emit: smorfs
    tuple val(meta), path("*/*.all_orfs.faa.gz")    , emit: all_orfs
    tuple val(meta), path("*/*.prediction.gz")      , emit: amp_prediction
    tuple val(meta), path("*/*.md")                 , emit: readme_file
    tuple val(meta), path("*/*_log.txt")            , emit: log_file
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macrel contigs \\
        $args \\
        --fasta $fasta \\
        --output ${prefix}/ \\
        --tag ${prefix} \\
        --log-file ${prefix}/${prefix}_log.txt \\
        --threads $task.cpus

    gzip --no-name ${prefix}/*.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macrel: \$(echo \$(macrel --version | sed 's/macrel //g'))
    END_VERSIONS
    """
}
