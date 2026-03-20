process WHATSHAP_PHASE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d837709891c2d98fc0956f6fd0dba18b0f67d96c4db74ccbae7db98fd00afe42/data'
        : 'community.wave.seqera.io/library/whatshap:2.8--7fe530bc624a3e5a' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(bam), path(bai)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('whatshap'), eval("whatshap --version"), emit: versions_whatshap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("${vcf}" == "${prefix}.vcf" || "${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    whatshap \\
        phase \\
        --output ${prefix}.vcf \\
        --reference ${fasta} \\
        ${args} \\
        ${vcf} \\
        ${bam}

    bgzip \\
        -@ ${task.cpus} \\
        ${prefix}.vcf

    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1 | sed 's/whatshap //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("${vcf}" == "${prefix}.vcf" || "${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
