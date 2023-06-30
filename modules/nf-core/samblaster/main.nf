process SAMBLASTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samblaster=0.1.26 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:cee56b506ceb753d4bbef7e05b81e1bfc25d937f-0' :
        'biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:cee56b506ceb753d4bbef7e05b81e1bfc25d937f-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$bam" == "${prefix}.bam" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools view -h $args2 $bam | \\
    samblaster $args | \\
    samtools view $args3 -Sb - >${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samblaster: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
