process PARACLU {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::paraclu=10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1' :
        'biocontainers/paraclu:10--h9a82719_1' }"

    input:
    tuple val(meta), path(bed)
    val(min_cluster)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '10' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    awk -F "\t" '{print\$1"\t"\$6"\t"\$2"\t"\$5}' < $bed > ${bed}_4P
    sort -k1,1 -k3n ${bed}_4P > ${bed}_4Ps
    paraclu $min_cluster ${bed}_4Ps > ${prefix}.clustered
    paraclu-cut  ${prefix}.clustered >  ${prefix}.clustered.simplified
    awk -F '\t' '{print \$1"\t"\$3"\t"\$4"\t"\$1":"\$3".."\$4","\$2"\t"\$6"\t"\$2}' ${prefix}.clustered.simplified >  ${prefix}.clustered.simplified.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paraclu: $VERSION
    END_VERSIONS
    """
}
