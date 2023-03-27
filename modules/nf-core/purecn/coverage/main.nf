process PURECN_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-purecn=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam)

    path intervals

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
   Rscript Coverage.R \
        --out-dir /purecn/coverage/${meta.id}/ \
        --bam ${bam} \
        --intervals ${intervals}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript /usr/local/lib/R/library/PureCN/extdata/PureCN.R --version)
    END_VERSIONS
    """
}
