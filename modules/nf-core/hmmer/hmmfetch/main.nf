process HMMER_HMMFETCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b171806a7b780537073c1e9ffe50885ba7917d8072b5941a125698df1253c4f7/data'
:         'community.wave.seqera.io/library/hmmer:3.4--ab3efa7c16f46920' }"

    // The module can be called with either a key, a file containing keys or neither.
    // In the latter case, the hmm database will be indexed and an index but no output
    // hmm will be produced.
    input:
    tuple val(meta), path(hmm)
    val   key
    path  keyfile
    path  index         // Only used to stage the index from a previous run

    output:
    tuple val(meta), path("*.hmm"), emit: hmm,   optional: true
    tuple val(meta), path("*.ssi"), emit: index, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def keyarg  = key ?: ''
    def kfopt   = keyfile ? '-f' : ''
    def index_opt   = ! key && ! keyfile ? '--index' : ''
    def outfile = ! key && ! keyfile ? '' : "> ${prefix}.hmm"

    """
    hmmfetch \\
        $kfopt \\
        $index_opt \\
        $args \\
        $hmm \\
        $keyarg \\
        $keyfile \\
        $outfile
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.hmm
    """
}
