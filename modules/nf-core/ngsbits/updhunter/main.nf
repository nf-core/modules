process NGSBITS_UPDHUNTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2be56a07ac1d5a447a10fd061be4d6144620bec00bac834f58c2bdef0330147f/data'
        : 'community.wave.seqera.io/library/ngs-bits:2025_09--f6ea3a4494373ed6'}"

    input:
    tuple val(meta), path(vcf), path(bed)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.igv"), emit: igv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_informative = args.contains('-out_informative') ? '' : "-out_informative ${prefix}.igv"
    def exclude = bed ? "-exclude ${bed}" : ''

    """
    UpdHunter \\
        -in ${vcf} \\
        ${exclude} \\
        ${out_informative} \\
        ${args} \\
        -out ${prefix}.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(UpdHunter --version 2>&1) | sed 's/UpdHunter //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.tsv
    touch ${prefix}.igv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(UpdHunter --version 2>&1) | sed 's/UpdHunter //' )
    END_VERSIONS
    """
}
