process SENTIEON_WGSMETRICS {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80ccb05eb4f1a193a3bd99c4da90f55f74ea6556c25f154e53e1ff5a6caa372d/data' :
        'community.wave.seqera.io/library/sentieon:202503--5e378058d837c58c' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(intervals_list)

    output:
    tuple val(meta), path('*.txt'), emit: wgs_metrics
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def interval = intervals_list  ? "--interval ${intervals_list}" : ""
    def input    = bam.sort().collect{"-i $it"}.join(' ')
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""
    """
    $sentieonLicense

    sentieon \\
        driver \\
        -t $task.cpus \\
        -r $fasta \\
        $input \\
        $interval \\
        $args \\
        --algo WgsMetricsAlgo \\
        ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
