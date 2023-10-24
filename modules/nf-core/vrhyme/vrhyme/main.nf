process VRHYME_VRHYME {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::vrhyme=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1':
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("vRhyme_best_bins_fasta/")                , emit: bins
    tuple val(meta), path("**/vRhyme_best_bins.*.membership.tsv")   , emit: membership
    tuple val(meta), path("**/vRhyme_best_bins.*.summary.tsv")      , emit: summary
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vRhyme \\
        -i $fasta \\
        -r $reads \\
        -o $prefix \\
        -t $task.cpus \\
        $args

    mv $prefix/vRhyme_best_bins_fasta/ vRhyme_best_bins_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    touch $prefix/vRhyme_best_bins.19.membership.tsv
    touch $prefix/vRhyme_best_bins.19.summary.tsv
    mkdir -p vRhyme_best_bins_fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_1.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_10.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_11.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_12.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_13.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_14.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_2.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_3.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_4.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_5.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_6.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_7.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_8.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_9.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' ))
    END_VERSIONS
    """
}
