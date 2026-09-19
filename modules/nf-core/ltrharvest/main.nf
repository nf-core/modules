process LTRHARVEST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ltr_harvest_parallel:1.1--hdfd78af_0':
        'quay.io/biocontainers/ltr_harvest_parallel:1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gff3") , emit: gff3
    tuple val(meta), path("*.scn")  , emit: scn
    tuple val("${task.process}"), val("LTR_HARVEST_parallel"), eval("LTR_HARVEST_parallel -h 2>&1 | sed -n 's/Version: v//p'"), emit: versions_ltr_harvest_parallel, topic: versions
    tuple val("${task.process}"), val("genometools"), eval("gt --version | sed -n 's/gt (GenomeTools) //p'"), emit: versions_genometools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    LTR_HARVEST_parallel \\
        -seq ${fasta} \\
        ${args} \\
        -threads ${task.cpus}

    mv "${fasta}.harvest.combine.gff3" \\
        "${prefix}.gff3"

    mv "${fasta}.harvest.combine.scn" \\
        "${prefix}.scn"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gff3"
    touch "${prefix}.scn"
    """
}
