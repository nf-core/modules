process GATK4_COMPOSESTRTABLEFILE {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c5ad1d2671d349f609383068e987b75faf2003404df3fff475ec80ba9a69c68e/data'
        : 'community.wave.seqera.io/library/gatk4-lite:4.7.0.0--ccccf4acf9b6eba5'}"

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
