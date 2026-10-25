process TCOFFEE_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e85c972da17c3e7d5e8b283cf61b106111180e3fc02d4b04f0933106b4ca50b3/data'
        : 'community.wave.seqera.io/library/t-coffee_pigz:aa6b4f1ff4adb194'}"


    input:
    tuple val(meta) , path(aln)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("*.{aln,aln.gz}")          , emit: alignment
    tuple val(meta), path("*.{score_html,sp_ascii}") , emit: eval, optional: true
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions
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
