process GATK4_INTERVALLISTTOOLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data':
        'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01' }"

    input:
    tuple val(meta), path(intervals)

    output:
    tuple val(meta), path("*_split/*/*.interval_list"), emit: interval_list
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK IntervalListTools] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    mkdir ${prefix}_split

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        IntervalListTools \\
        --INPUT ${intervals} \\
        --OUTPUT ${prefix}_split \\
        --TMP_DIR . \\
        ${args}

    python3 <<CODE
    import glob, os
    # The following python code snippet rename the output files into different name to avoid overwriting or name conflict
    intervals = sorted(glob.glob("*_split/*/*.interval_list"))
    for i, interval in enumerate(intervals):
        (directory, filename) = os.path.split(interval)
        newName = os.path.join(directory, str(i + 1) + filename)
        os.rename(interval, newName)
    CODE
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_split/temp_0001_of_6
    mkdir -p ${prefix}_split/temp_0002_of_6
    mkdir -p ${prefix}_split/temp_0003_of_6
    mkdir -p ${prefix}_split/temp_0004_of_6
    touch ${prefix}_split/temp_0001_of_6/1scattered.interval_list
    touch ${prefix}_split/temp_0002_of_6/2scattered.interval_list
    touch ${prefix}_split/temp_0003_of_6/3scattered.interval_list
    touch ${prefix}_split/temp_0004_of_6/4scattered.interval_list
    """
}
