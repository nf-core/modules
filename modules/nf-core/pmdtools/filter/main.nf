process PMDTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pmdtools=0.60" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmdtools:0.60--hdfd78af_5' :
        'quay.io/biocontainers/pmdtools:0.60--hdfd78af_5' }"

    input:
    tuple val(meta), path(bam), path (bai)
    val(threshold)
    path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def split_cpus = Math.floor(task.cpus/2)
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "[pmdtools/filter] Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    //threshold and header flags activate filtering function of pmdtools
    """
    samtools \\
        calmd \\
        $bam \\
        $reference \\
        $args \\
        -@ ${split_cpus} \\
    | pmdtools \\
        --threshold $threshold \\
        --header \\
        $args2 \\
    | samtools \\
        view \\
        $args3 \\
        -Sb \\
        - \\
        -@ ${split_cpus} \\
        -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmdtools: \$( pmdtools --version | cut -f2 -d ' ' | sed 's/v//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
