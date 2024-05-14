process BISCUIT_BSCONV {
    tag "$meta.id"
    label 'process_long'


    conda "bioconda::biscuit=1.1.0.20220707"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.1.0.20220707--he272189_1':
        'biocontainers/biscuit:1.1.0.20220707--he272189_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bsconv_bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/\\.bis.amb\$//'`

    biscuit bsconv \\
        $args \\
        \$INDEX \\
        $bam \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}
