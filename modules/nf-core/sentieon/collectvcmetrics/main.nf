process SENTIEON_COLLECTVCMETRICS {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a64461f38d76bebea8e21441079e76e663e1168b0c59dafee6ee58440ad8c8ac/data' :
        'community.wave.seqera.io/library/sentieon:202308.03--59589f002351c221' }"

    input:
    tuple val(meta) , path(vcf), path(tbi)
    tuple val(meta2), path(dbsnp), path(dbsnp_tbi)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(interval)

    output:
    tuple val(meta), path("*.variant_calling_detail_metrics") , emit: metrics
    tuple val(meta), path("*.variant_calling_summary_metrics"), emit: summary
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def args         = task.ext.args  ?: ''
    def args2        = task.ext.args2 ?: ''
    def interval_cmd = interval ? "--interval $interval" : ""
    """
    sentieon \\
        driver \\
        -t $task.cpus \\
        -r $fasta \\
        $interval_cmd \\
        $args \\
        --algo CollectVCMetrics \\
        -v ${vcf} \\
        -d ${dbsnp} \\
        $args2 \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_calling_detail_metrics
    touch ${prefix}.variant_calling_summary_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
