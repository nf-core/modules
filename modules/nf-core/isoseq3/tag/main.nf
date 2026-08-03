process ISOSEQ3_TAG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:4.0.0--h9ee0642_0':
        'quay.io/biocontainers/isoseq:4.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    val design

    output:
    tuple val(meta), path("*.flt.bam")                  , emit: bam
    tuple val(meta), path("*.flt.bam.pbi")              , emit: pbi
    tuple val("${task.process}"), val('isoseq3'), eval("isoseq tag --version | head -n 1 | sed 's/isoseq tag //g' | sed 's/ (.*//g'"), emit: versions_isoseq3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_design = ~/^(?:\d{1,2}[UBGX]-)+T$|^(?:\d{1,2}[UBGX]-)+T(?:-\d{1,2}[UBGX])+$|^T(?:-\d{1,2}[UBGX])+$/
    if ( !(design ==~ valid_design) )  { error "Invalid UMI/barcode design. Check https://isoseq.how/umi/umi-barcode-design.html for how to specify the design" }
    """
    isoseq \\
        tag \\
        -j $task.cpus \\
        --design $design \\
        ${prefix}.5p--3p.bam \\
        ${prefix}.flt.bam \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.flt.bam
    touch ${prefix}.flt.bam.pbi
    """
}
