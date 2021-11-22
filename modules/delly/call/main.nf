process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::delly=0.8.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/delly:0.8.7--he03298f_1"
    } else {
        container "quay.io/biocontainers/delly:0.8.7--he03298f_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val(meta), path("*.csi"), emit: csi
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    delly \\
        call \\
        $args \\
        -o ${prefix}.bcf \\
        -g  $fasta \\
        $bam \\

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
