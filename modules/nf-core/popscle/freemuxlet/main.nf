process POPSCLE_FREEMUXLET {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popscle:0.1beta--h2c78cec_0' :
        'biocontainers/popscle:0.1beta--h2c78cec_0' }"

    input:
    tuple val(meta), path(plp), val(n_sample)

    output:
    tuple val(meta), path('*.clust1.samples.gz')  , emit: result
    tuple val(meta), path('*.clust1.vcf.gz')      , emit: vcf
    tuple val(meta), path('*.lmix')               , emit: lmix
    tuple val(meta), path('*.clust0.samples.gz')  , emit: singlet_result   , optional: true
    tuple val(meta), path('*.clust0.vcf.gz')      , emit: singlet_vcf      , optional: true
    path 'versions.yml'                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    popscle freemuxlet \\
        --plp ${plp}/$prefix \\
        --out $prefix \\
        --nsample $n_sample \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    touch ${prefix}.clust1.samples.gz
    touch ${prefix}.clust1.vcf.gz
    touch ${prefix}.lmix

    if [[ "$args" == *"--aux-files"* ]]; then
        touch ${prefix}.clust0.samples.gz
        touch ${prefix}.clust0.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle: $VERSION
    END_VERSIONS
    """
}
