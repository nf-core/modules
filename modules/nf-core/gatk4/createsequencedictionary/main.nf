process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${fasta}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/06957b3d467331a29a6d770308f17004738469a21c0aedf893d5ea0312accaea/data'
        : 'community.wave.seqera.io/library/gatk4-main:4.7.0.0--b8e8676903b345d8'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.dict'), emit: dict
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def avail_mem = 6144
    if (!task.memory) {
        log.info('[GATK CreateSequenceDictionary] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CreateSequenceDictionary \\
        --REFERENCE ${fasta} \\
        --URI ${fasta} \\
        --TMP_DIR . \\
        ${args}
    """

    stub:
    """
    touch ${fasta.baseName}.dict
    """
}
