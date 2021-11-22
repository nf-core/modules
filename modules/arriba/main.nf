process ARRIBA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::arriba=2.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/arriba:2.1.0--h3198e80_1"
    } else {
        container "quay.io/biocontainers/arriba:2.1.0--h3198e80_1"
    }

    input:
    tuple val(meta), path(bam)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*.fusions.tsv")          , emit: fusions
    tuple val(meta), path("*.fusions.discarded.tsv"), emit: fusions_fail
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def blacklist = (options.args.contains('-b')) ? '' : '-f blacklist'
    """
    arriba \\
        -x $bam \\
        -a $fasta \\
        -g $gtf \\
        -o ${prefix}.fusions.tsv \\
        -O ${prefix}.fusions.discarded.tsv \\
        $blacklist \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
