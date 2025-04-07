process RIBOTRICER_PREPAREORFS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotricer:1.3.3--pyhdfd78af_0':
        'biocontainers/ribotricer:1.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("*_candidate_orfs.tsv"), emit: candidate_orfs
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ribotricer prepare-orfs \\
        --gtf $gtf \\
        --fasta $fasta \\
        --prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotricer: \$(ribotricer --version | grep ribotricer |& sed '1!d ; s/ribotricer, version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_candidate_orfs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotricer: \$(ribotricer --version | grep ribotricer |& sed '1!d ; s/ribotricer, version //')
    END_VERSIONS
    """
}
