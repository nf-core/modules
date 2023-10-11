process HMMER_HMMFETCH {
    tag "$hmm"
    label 'process_single'

    conda "bioconda::hmmer=3.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h87f3376_2':
        'biocontainers/hmmer:3.3.2--h87f3376_2' }"

    // The module can be called with either a key, a file containing keys or neither.
    // In the latter case, the hmm database will be indexed and an index but no output
    // hmm will be produced.
    input:
    path  hmm
    val   key
    path  keyfile
    path  index         // Only used to stage the index from a previous run

    output:
    path "selection.hmm", emit: hmm,   optional: true
    path "*.ssi"        , emit: index, optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def keyarg  = key ?: ''
    def kfopt   = keyfile ? '-f' : ''
    def index   = ! key && ! keyfile ? '--index' : ''
    def outfile = ! key && ! keyfile ? '' : '> selection.hmm'

    """
    hmmfetch \\
        $kfopt \\
        $index \\
        $args \\
        $hmm \\
        $keyarg \\
        $keyfile \\
        $outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch selection.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
