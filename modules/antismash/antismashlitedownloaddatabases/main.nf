process ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES {
    tag '$database_dir'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::antismash-lite=6.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash-lite:6.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash-lite:6.0.1--pyhdfd78af_0' }"

    input:

    output:
    path("*/clusterblast")   , emit: database_clusterblast
    path("*/clustercompare") , emit: database_clustercompare
    path("*/pfam")           , emit: database_pfam
    path("*/resfam")         , emit: database_resfam
    path("*/tigrfam")        , emit: database_tigrfam
    path("*/README")         , emit: database_readme
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    download-antismash-databases \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
