process BISMARK_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe2a9d58209a38df5c99615a41d9ed6e8e546380d04c176e076e107010819a72/data' :
        'community.wave.seqera.io/library/bismark:0.25.0--95ba99b483e2eaf9' }"

    input:
    val(bam)
    path(align_report)
    path(dedup_report)
    path(splitting_report)
    path(mbias)

    output:
    path("*report.{html,txt}"), emit: summary
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark2summary ${bam.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bam.baseName()}"
    """
    touch ${prefix}.report.txt
    touch ${prefix}.report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
