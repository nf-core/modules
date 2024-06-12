process METAMAPS_MAPDIRECTLY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metamaps:0.1.633d2e0--h21ec9f0_0':
        'biocontainers/metamaps:0.1.633d2e0--h21ec9f0_0' }"

    input:
    tuple val(meta), path(reads)
    path database

    output:
    tuple val(meta), path("*classification_res")                          , emit: classification_res
    tuple val(meta), path("*classification_res.meta")                     , emit: meta_file
    tuple val(meta), path("*classification_res.meta.unmappedReadsLengths"), emit: meta_unmappedreadsLengths
    tuple val(meta), path("*classification_res.parameters")               , emit: para_file
    path "versions.yml"                                                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db=`find -L ${database} -name "DB.fa"`
    metamaps \\
        mapDirectly \\
        $args \\
        --all \\
        --reference \$db \\
        --threads $task.cpus \\
        --query $reads \\
        --output ${prefix}.classification_res

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps: \$(metamaps | sed -n 2p | sed 's/^.*MetaMaps v //')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_classification_res
    touch ${prefix}_classification_res.meta
    touch ${prefix}_classification_res.meta.unmappedReadsLengths
    touch ${prefix}_classification_res.parameters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps: \$(metamaps | sed -n 2p | sed 's/^.*MetaMaps v //')
    END_VERSIONS
    """
}
