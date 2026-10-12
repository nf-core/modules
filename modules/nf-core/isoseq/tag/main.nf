process ISOSEQ_TAG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:26.2.0--h9ee0642_0':
        'quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    val design

    output:
    tuple val(meta), path("*.flt.bam")    , emit: bam
    tuple val(meta), path("*.flt.bam.pbi"), emit: pbi
    tuple val("${task.process}"), val('isoseq'), eval("isoseq tag --version | sed -n '1s/isoseq tag \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_isoseq

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
        $bam \\
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
