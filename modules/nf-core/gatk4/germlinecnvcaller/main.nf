process GATK4_GERMLINECNVCALLER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(tsv), path(intervals), path(ploidy), path(model)

    output:
    tuple val(meta), path("*-cnv-model/*-calls"), emit: cohortcalls, optional: true
    tuple val(meta), path("*-cnv-model/*-model"), emit: cohortmodel, optional: true
    tuple val(meta), path("*-cnv-calls/*-calls"), emit: casecalls, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals ${intervals}" : ""
    def ploidy_command = ploidy ? "--contig-ploidy-calls ${ploidy}" : ""
    def model_command = model ? "--model ${model}" : ""
    def input_list = tsv.collect { tsv_ -> "--input ${tsv_}" }.join(' ')
    def output_command = model ? "--output ${prefix}-cnv-calls" : "--output ${prefix}-cnv-model"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GermlineCNVCaller \\
        ${input_list} \\
        ${ploidy_command} \\
        ${output_command} \\
        --output-prefix ${prefix} \\
        ${args} \\
        ${intervals_command} \\
        ${model_command}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}-cnv-calls/${prefix}-calls
    mkdir -p ${prefix}-cnv-model/${prefix}-model
    mkdir -p ${prefix}-cnv-model/${prefix}-calls
    """
}
