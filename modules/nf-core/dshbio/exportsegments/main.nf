process DSHBIO_EXPORTSEGMENTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dsh-bio:3.0--hdfd78af_0' :
        'quay.io/biocontainers/dsh-bio:3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.fa.gz"), emit: fasta
    tuple val("${task.process}"), val('dsh-bio'), eval("dsh-bio --version | sed '1!d;s/dsh-bio-tools //'"), emit: versions_dshbio, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dsh-bio \\
        export-segments \\
        $args \\
        -i $gfa \\
        -o ${prefix}.fa.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.fa.gz
    """
}
