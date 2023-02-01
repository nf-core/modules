
process VGAN_HAPLOCART {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vgan=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vgan:1.0.1--h9ee0642_0':
        'quay.io/biocontainers/vgan:1.0.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(hc_files), path(reads), path(reads2)

    output:
    tuple val(meta), path("*[!posterior].txt") , emit: txt   
    tuple val(meta), path("*.posterior.txt)    , emit: posterior
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1'
    def reads_args = meta.single_end ? "-fq1 ${reads}" : "-fq1 ${reads} -fq2 ${reads2}"

    """
    vgan haplocart \\
        $args \\
        -t $task.cpus \\
        $reads_args \\
        -o ${prefix}.txt \\
        --hc-files ${hc_files} \\
        -pf ${prefix}.posterior.txt

    echo $VERSION >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.posterior.txt

    echo $VERSION >> versions.yml
    """

}
