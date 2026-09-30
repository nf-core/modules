process GATK4_COMPOSESTRTABLEFILE {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd65ba227372b7f3f62231d105427edcf5def4c4700689f77ce9536c954a22f6/data'
        : 'community.wave.seqera.io/library/gatk4-main:4.7.0.0--6f748daadc3eeb06'}"

    input:
    path fasta
    path fasta_fai
    path dict

    output:
    path "*.zip", emit: str_table
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def avail_mem = 6144
    if (!task.memory) {
        log.info('[GATK ComposeSTRTableFile] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ComposeSTRTableFile \\
        --reference ${fasta} \\
        --output ${fasta.baseName}.zip \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    """
    touch ${fasta.baseName}.zip
    """
}
