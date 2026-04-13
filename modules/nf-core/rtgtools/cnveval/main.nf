process RTGTOOLS_CNVEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(query_vcf), path(query_vcf_tbi), path(truth_vcf), path(truth_vcf_tbi), path(evaluation_regions_bed)

    output:
    tuple val(meta), path("*.baseline.bed.gz")    , emit: baseline_bed
    tuple val(meta), path("*.calls.bed.gz")       , emit: calls_bed
    tuple val(meta), path("*.weighted_roc.tsv.gz"), emit: weighted_roc
    tuple val(meta), path("*.summary.txt")        , emit: summary
    tuple val("${task.process}"), val('rtgtools'), eval("rtg version | sed 's/Product: RTG Tools //; q'"), topic: versions, emit: versions_rtgtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory.toGiga() + "G"

    """
    rtg RTG_MEM=$avail_mem cnveval \\
        ${args} \\
        --baseline=${truth_vcf} \\
        --calls=${query_vcf} \\
        --evaluation-regions=${evaluation_regions_bed} \\
        --output=${prefix} \\
        --threads=${task.cpus}

    mv ${prefix}/done ${prefix}/progress . 2>/dev/null || true
    for f in ${prefix}/*; do
        mv "\$f" "${prefix}.\${f##*/}"
    done
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo | gzip > ${prefix}.baseline.bed.gz
    echo | gzip > ${prefix}.calls.bed.gz
    echo | gzip > ${prefix}.weighted_roc.tsv.gz
    touch ${prefix}.summary.txt
    """
}
