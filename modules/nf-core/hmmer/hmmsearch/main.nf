process HMMER_HMMSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(hmmfile), path(seqdb), val(write_align), val(write_target), val(write_domain)

    output:
    tuple val(meta), path('*.txt.gz')   , emit: output
    tuple val(meta), path('*.sto.gz')   , emit: alignments    , optional: true
    tuple val(meta), path('*.tbl.gz')   , emit: target_summary, optional: true
    tuple val(meta), path('*.domtbl.gz'), emit: domain_summary, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output         = "${prefix}.txt"
    def alignment      = write_align  ? "-A ${prefix}.sto" : ''
    def target_summary = write_target ? "--tblout ${prefix}.tbl" : ''
    def domain_summary = write_domain ? "--domtblout ${prefix}.domtbl" : ''
    def to_gzip = ["${prefix}.txt"]
    if (write_align)  to_gzip << "${prefix}.sto"
    if (write_target) to_gzip << "${prefix}.tbl"
    if (write_domain) to_gzip << "${prefix}.domtbl"
    """
    hmmsearch \\
        $args \\
        --cpu $task.cpus \\
        -o $output \\
        $alignment \\
        $target_summary \\
        $domain_summary \\
        $hmmfile \\
        $seqdb

    ${to_gzip ? "gzip --no-name ${to_gzip.join(' ')}" : ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def to_gzip = ["${prefix}.txt"]
    if (write_align)  to_gzip << "${prefix}.sto"
    if (write_target) to_gzip << "${prefix}.tbl"
    if (write_domain) to_gzip << "${prefix}.domtbl"
    """
    ${to_gzip ? "touch ${to_gzip.join(' ')}" : ''}

    ${to_gzip ? "gzip --no-name ${to_gzip.join(' ')}" : ''}
    """
}
