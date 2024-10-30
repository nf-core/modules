process SCDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-scds:1.18.0--r43hdfd78af_0':
        'biocontainers/bioconductor-scds:1.18.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.csv"), emit: predictions
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'scds.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds
    touch ${prefix}.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scds: \$(Rscript -e "library(scds); cat(as.character(packageVersion('scds')))")
    END_VERSIONS
    """
}
