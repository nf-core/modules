process VARSCAN_FPFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed57a091507c62e990bbd08d532281d161d99f060316e0a991791f167d7b1daf/data':
        'community.wave.seqera.io/library/htslib_varscan:24b3b3db2ca78de8' }"

    input:
    tuple val(meta), path(vcf), path(rc)

    output:
    tuple val(meta), path("*.pass.vcf.gz"), emit: pass_vcf
    tuple val(meta), path("*.fail.vcf.gz"), emit: fail_vcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def compressed = vcf.name.endsWith('.gz') ? true : false
    """
    ${compressed ? "bgzip -df $vcf" : ""}

    varscan fpfilter \\
        ${compressed ? "${vcf.baseName}" : "$vcf"} \\
        $rc \\
        --output-file ${prefix}.pass.vcf \\
        --filtered-file ${prefix}.fail.vcf \\
        $args

    bgzip ${prefix}.pass.vcf
    bgzip ${prefix}.fail.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    echo "" | gzip > ${prefix}.pass.vcf.gz
    echo "" | gzip > ${prefix}.fail.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
