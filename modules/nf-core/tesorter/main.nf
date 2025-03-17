process TESORTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tesorter:1.4.6--pyhdfd78af_0':
        'biocontainers/tesorter:1.4.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db_hmm

    output:
    tuple val(meta), path("*.domtbl")   , emit: domtbl
    tuple val(meta), path("*.dom.faa")  , emit: dom_faa
    tuple val(meta), path("*.dom.tsv")  , emit: dom_tsv
    tuple val(meta), path("*.dom.gff3") , emit: dom_gff3
    tuple val(meta), path("*.cls.tsv")  , emit: cls_tsv , optional: true
    tuple val(meta), path("*.cls.lib")  , emit: cls_lib , optional: true
    tuple val(meta), path("*.cls.pep")  , emit: cls_pep , optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def db_hmm_arg  = db_hmm ? "--db-hmm $db_hmm" : ''
    """
    TEsorter \\
        $db_hmm_arg \\
        --processors $task.cpus \\
        -pre $prefix \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def touch_opts  = ( ! args.contains('-genome') ) ? 'Yes' : 'No'
    """
    touch ${prefix}.domtbl
    touch ${prefix}.dom.faa
    touch ${prefix}.dom.tsv
    touch ${prefix}.dom.gff3

    if [ "$touch_opts" = "Yes" ]; then
        touch ${prefix}.cls.tsv
        touch ${prefix}.cls.lib
        touch ${prefix}.cls.pep
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
    END_VERSIONS
    """
}
