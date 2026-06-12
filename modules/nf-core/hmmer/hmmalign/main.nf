process HMMER_HMMALIGN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(fasta)
    path hmm

    output:
    tuple val(meta), path("*.sto.gz"), emit: sto
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hmmalign \\
        $args \\
        $hmm \\
        $fasta | gzip -c > ${prefix}.sto.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.sto.gz
    """
}
