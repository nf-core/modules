process HMMER_HMMBUILD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(alignment)
    path mxfile

    output:
    tuple val(meta), path("*.hmm.gz"), emit: hmm
    path "*.hmmbuild.txt",             emit: hmmbuildout
    tuple val("${task.process}"), val("hmmer"), eval("hmmbuild -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER //'"), topic: versions, emit: versions_hmmer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def mxfileopt = mxfile ? "--mxfile ${mxfile}" : ""

    """
    hmmbuild \\
        $args \\
        --cpu $task.cpus \\
        -n ${prefix}  \\
        -o ${prefix}.hmmbuild.txt \\
        ${mxfileopt} \\
        ${prefix}.hmm \\
        $alignment

    gzip ${prefix}.hmm
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.hmm.gz
    touch ${prefix}.hmmbuild.txt
    """
}
