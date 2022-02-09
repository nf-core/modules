process HMMER_HMMSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1' :
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    tuple val(meta), path(hmmfile), path(seqdb), val(write_align), val(write_target), val(write_domain)

    output:
    tuple val(meta), path(output)    , emit: output
    tuple val(meta), path('*.sto')   , emit: alignments    , optional: true
    tuple val(meta), path('*.tbl')   , emit: target_summary, optional: true
    tuple val(meta), path('*.domtbl'), emit: domain_summary, optional: true
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.txt"
    alignment = write_align ? "-A ${prefix}.sto" : ''
    target_summary = write_target ? "--tblout ${prefix}.tbl" : ''
    domain_summary = write_domain ? "--domtblout ${prefix}.domtbl" :  ''
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
