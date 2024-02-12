
process UMITOOLS_PREPAREFORRSEM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::umi_tools=1.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.4--py38hbff2b2d_1' :
        'biocontainers/umi_tools:1.1.4--py38hbff2b2d_1' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'prepare-for-rsem.py'

    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
