process ALLELECOUNTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.3.0--h41abebc_0' :
        'biocontainers/cancerit-allelecount:4.3.0--h41abebc_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path loci
    path fasta

    output:
    tuple val(meta), path("*.alleleCount"), emit: allelecount
    tuple val("${task.process}"), val('alleleCounter'), eval('alleleCounter --version'), emit: versions_allelecounter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_options = fasta ? "-r $fasta": ""
    """
    alleleCounter \\
        ${args} \\
        -l ${loci} \\
        -b ${input} \\
        ${reference_options} \\
        -o ${prefix}.alleleCount
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alleleCount
    """
}
