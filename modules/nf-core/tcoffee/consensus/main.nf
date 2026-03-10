process TCOFFEE_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5c1c07cc0099dacea172348bc78f9a9baab592ce3ece89873703b9e963d269/data':
        'community.wave.seqera.io/library/t-coffee_pigz:c98a6c87c62d9df6' }"


    input:
    tuple val(meta) , path(aln)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.{aln,aln.gz}")          , emit: alignment
    tuple val(meta), path("*.{score_html,sp_ascii}") , emit: eval, optional: true
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | awk \'{gsub("Version_", ""); print \\$3}\''), emit: versions_tcoffee, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tree_args = tree ? "-usetree $tree" : ""
    def outfile = compress ? "stdout" : "${prefix}.aln"
    def write_output = compress ? " | pigz -cp ${task.cpus} > ${prefix}.aln.gz" : ""
    """
    export TEMP='./'
    t_coffee -aln ${aln} \
        $tree_args \
        $args \
        -thread ${task.cpus} \
        -outfile $outfile \
        $write_output
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TEMP='./'
    touch ${prefix}.aln${compress ? '.gz':''}
    """
}
