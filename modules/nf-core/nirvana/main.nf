process NIRVANA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/nirvana:3.18.1--910f092f78f85c70'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(cache), val(cache_prefix)
    tuple val(meta4), path(supplementary_annotations)

    output:
    tuple val(meta), path("*.json.gz"), emit: json
    tuple val(meta), path("*.json.gz.jsi"), emit: jsi
    tuple val("${task.process}"), val('nirvana'), eval("Nirvana -v 2>&1 | awk '{print \$2}' | cut -d'-' -f1"), topic: versions, emit: versions_nirvana

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-c ${cache}/${cache_prefix}" : ""
    def sa_command = supplementary_annotations ? "--sd ${supplementary_annotations}" : ""
    """
    Nirvana \\
        -i ${vcf} \\
        -r ${reference} \\
        ${cache_command} \\
        ${sa_command} \\
        ${args} \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.json.gz
    touch ${prefix}.json.gz.jsi
    """
}
