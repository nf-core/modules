process PICARD_COLLECTHSMETRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(bam), path(bai), path(bait_intervals, stageAs: "baits/*"), path(target_intervals, stageAs: 'targets/*')
    tuple val(meta2), path(ref)
    tuple val(meta3), path(ref_fai)
    tuple val(meta4), path(ref_dict)
    tuple val(meta5), path(ref_gzi)

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val("${task.process}"), val('picard'), eval("picard CollectHsMetrics --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = ref ? "--REFERENCE_SEQUENCE ${ref}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    def bait_interval_list = bait_intervals
    def bait_intervallist_cmd = ""
    if (bait_intervals =~ /.(bed|bed.gz)$/) {
        bait_interval_list = bait_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        bait_intervallist_cmd = "picard -Xmx${avail_mem}M BedToIntervalList --INPUT ${bait_intervals} --OUTPUT ${bait_interval_list} --SEQUENCE_DICTIONARY ${ref_dict} --TMP_DIR ."
    }

    def target_interval_list = target_intervals
    def target_intervallist_cmd = ""
    if (target_intervals =~ /.(bed|bed.gz)$/) {
        target_interval_list = target_intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        target_intervallist_cmd = "picard -Xmx${avail_mem}M BedToIntervalList --INPUT ${target_intervals} --OUTPUT ${target_interval_list} --SEQUENCE_DICTIONARY ${ref_dict} --TMP_DIR ."
    }
    """
    ${bait_intervallist_cmd}
    ${target_intervallist_cmd}

    picard \\
        -Xmx${avail_mem}M \\
        CollectHsMetrics \\
        ${args} \\
        ${reference} \\
        --BAIT_INTERVALS ${bait_interval_list} \\
        --TARGET_INTERVALS ${target_interval_list} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.CollectHsMetrics.coverage_metrics
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectHsMetrics.coverage_metrics
    """
}
