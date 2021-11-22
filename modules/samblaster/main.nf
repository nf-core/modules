process SAMBLASTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samblaster=0.1.26 bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:ba4a02b56f3e524a6e006bcd99fe8cc1d7fe09eb-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:ba4a02b56f3e524a6e006bcd99fe8cc1d7fe09eb-0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if( "$bam" == "${prefix}.bam" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    samtools view -h $task.ext.args2 $bam | \\
    samblaster $args | \\
    samtools view $task.ext.args3 -Sb - >${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
