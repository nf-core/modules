process GATK4_POSTPROCESSGERMLINECNVCALLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(calls), path(model), path(ploidy)

    output:
    tuple val(meta), path("*_genotyped_intervals.vcf.gz"), emit: intervals, optional: true
    tuple val(meta), path("*_genotyped_segments.vcf.gz"), emit: segments, optional: true
    tuple val(meta), path("*_denoised.vcf.gz"), emit: denoised, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def calls_command = calls ? calls.collect { call -> "--calls-shard-path ${call}" }.join(' ') : ""
    def model_command = model ? model.collect { mode -> "--model-shard-path ${mode}" }.join(' ') : ""
    def ploidy_command = ploidy ? "--contig-ploidy-calls ${ploidy}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PostProcessGermlineCnvCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PostprocessGermlineCNVCalls \\
        ${calls_command} \\
        ${model_command} \\
        ${ploidy_command} \\
        ${args} \\
        --output-genotyped-intervals ${prefix}_genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${prefix}_genotyped_segments.vcf.gz \\
        --output-denoised-copy-ratios ${prefix}_denoised.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_genotyped_intervals.vcf.gz
    echo "" | gzip > ${prefix}_genotyped_segments.vcf.gz
    echo "" | gzip > ${prefix}_denoised.vcf.gz
    """
}
