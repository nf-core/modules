process FGBIO_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87626ef674e2f19366ae6214575a114fe80ce598e796894820550731706a84be/data' :
        'community.wave.seqera.io/library/fgbio:2.4.0--913bad9d47ff8ddc' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    val(min_reads)
    val(min_baseq)
    val(max_base_error_rate)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_consensus_filtered"

    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio FilterConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --compression=0 \\
        FilterConsensusReads \\
        --input $bam \\
        --output ${prefix}.bam \\
        --ref ${fasta} \\
        --min-reads ${min_reads} \\
        --min-base-quality ${min_baseq} \\
        --max-base-error-rate ${max_base_error_rate} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_consensus_filtered"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
