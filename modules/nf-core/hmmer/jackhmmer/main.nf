process HMMER_JACKHMMER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b171806a7b780537073c1e9ffe50885ba7917d8072b5941a125698df1253c4f7/data'
:         'community.wave.seqera.io/library/hmmer:3.4--ab3efa7c16f46920' }"

    input:
    tuple val(meta), path(fasta), path(seqdb), val(write_align), val(write_target), val(write_domain)

    output:
    tuple val(meta), path('*.txt.gz')   , emit: output
    tuple val(meta), path('*.sto.gz')   , emit: alignments    , optional: true
    tuple val(meta), path('*.tbl.gz')   , emit: target_summary, optional: true
    tuple val(meta), path('*.domtbl.gz'), emit: domain_summary, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("jackhmmer -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    alignment      = write_align  ? "-A ${prefix}.sto" : ''
    target_summary = write_target ? "--tblout ${prefix}.tbl" : ''
    domain_summary = write_domain ? "--domtblout ${prefix}.domtbl" : ''
    """
    jackhmmer \\
        $args \\
        --cpu $task.cpus \\
        -o ${prefix}.txt \\
        $alignment \\
        $target_summary \\
        $domain_summary \\
        ${fasta} \\
        ${seqdb}

    gzip --no-name *.txt \\
        ${write_align ? '*.sto' : ''} \\
        ${write_target ? '*.tbl' : ''} \\
        ${write_domain ? '*.domtbl' : ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.txt"
    ${write_align ? "touch ${prefix}.sto" : ''} \\
    ${write_target ? "touch ${prefix}.tbl" : ''} \\
    ${write_domain ? "touch ${prefix}.domtbl" : ''}

    gzip --no-name *.txt \\
        ${write_align ? '*.sto' : ''} \\
        ${write_target ? '*.tbl' : ''} \\
        ${write_domain ? '*.domtbl' : ''}
    """
}
