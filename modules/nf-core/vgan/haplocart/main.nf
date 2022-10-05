def VERSION = '1.0.0'

process VGAN_HAPLOCART {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vgan=1.0.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vgan:1.0.0--h9ee0642_0':
        'quay.io/biocontainers/vgan' }"

    input:
    tuple val(meta), path(inputfile)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools bam2fq $inputfile | vgan \\
        haplocart \\
        $args \\
        -fq1 /dev/stdin \\
        -o ${prefix}.txt \\

    echo $VERSION >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    echo $VERSION >> versions.yml
    """

}
