process DYSGU_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1a8c4c95342498790fe752b702051dc40eb71114d4c0b36e844daad1fea7b593/data':
        'community.wave.seqera.io/library/dysgu:1.8.7--a06ec137d500dc83' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(sites)
    tuple val(meta5), path(bed)
    tuple val(meta6), path(search_bed)
    tuple val(meta7), path(exclude_bed)

    output:
    tuple val(meta), path('*.vcf.gz')       , emit: vcf
    tuple val(meta), path('*.vcf.gz.tbi')   , emit: tbi
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def sites_vcf = sites ? "--sites ${sites}" : ''
    def regions = bed ? "--regions ${bed}" : ''
    def search = search_bed ? "--search ${search_bed}" : ''
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ''

    """
    dysgu run \\
        ${args} \\
        ${sites_vcf} \\
        ${regions} \\
        ${search} \\
        ${exclude} \\
        --procs ${task.cpus} \\
        --overwrite \\
        $fasta \\
        . \\
        $input \\
        | bgzip ${args2} --threads ${task.cpus} --stdout > ${prefix}.vcf.gz && tabix ${args3} ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dysgu: \$(dysgu --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dysgu: \$(dysgu --version 2>&1)
    END_VERSIONS
    """
}
