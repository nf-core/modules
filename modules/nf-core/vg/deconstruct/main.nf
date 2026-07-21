process VG_DECONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.73.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.73.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(gfa)
    path(pb)
    path(gbwt)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('vg'), eval("vg 2>&1 | sed -n 's/.*version v\\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_vg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def snarls = pb ? "--snarls ${pb}" : ""
    def gbwt_arg = gbwt ? "--gbwt ${gbwt}" : ""
    """
    vg deconstruct \\
        --threads $task.cpus \\
        $args \\
        $snarls \\
        $gbwt_arg \\
        $gfa > ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
