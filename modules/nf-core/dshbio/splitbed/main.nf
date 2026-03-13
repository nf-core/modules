process DSHBIO_SPLITBED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dsh-bio:3.0--hdfd78af_0' :
        'biocontainers/dsh-bio:3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    tuple val("${task.process}"), val('dsh-bio'), eval("dsh-bio --version | head -n 1 | sed 's/dsh-bio-tools //'"), emit: versions_dshbio, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dsh-bio \\
        split-bed \\
        $args \\
        -p $prefix \\
        -s '.bed.gz' \\
        -i $bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}0.bed.gz
    echo | gzip > ${prefix}1.bed.gz
    """
}
