process SEMIBIN_SINGLEEASYBIN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::semibin=1.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/semibin%3A1.4.0--pyh7cba7a3_0':
        'quay.io/biocontainers/semibin%3A1.4.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta), path(bam)
    val min_len
    val ml_threshold
    val minfasta_kbs	

    output:
    tuple val(meta), path("*.csv")                        , emit: csv
    tuple val(meta), path("*.h5")                         , emit: model
    tuple val(meta), path("output_prerecluster_bins/*.fa"), emit: output_fasta
    tuple val(meta), path("output_recluster_bins/*.fa")   , emit: recluster_fasta
    tuple val(meta), path("*.tsv")                        , emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    SemiBin \\
	$args \\
        single_easy_bin \\
        -i $fasta \\
        -b $bam \\
        -o $prefix \\
        -t $task.cpus \\
        --min-len  $min_len \\
        --self-supervised  \\
        --ml-threshold $ml_threshold  \\
        --minfasta-kbs $minfasta_kbs       
    
    mv $prefix/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}": SemiBin: \$( SemiBin --version )
END_VERSIONS
"""
}
