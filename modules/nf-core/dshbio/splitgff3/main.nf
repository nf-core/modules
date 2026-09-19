process DSHBIO_SPLITGFF3 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dsh-bio:3.0--hdfd78af_0' :
        'quay.io/biocontainers/dsh-bio:3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gff3.gz"), emit: gff3
    tuple val("${task.process}"), val('dsh-bio'), eval("dsh-bio --version | sed '1!d;s/dsh-bio-tools //'"), emit: versions_dshbio, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dsh-bio \\
        split-gff3 \\
        $args \\
        -p $prefix \\
        -s '.gff3.gz' \\
        -i $gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.gff3.gz
    """
}
