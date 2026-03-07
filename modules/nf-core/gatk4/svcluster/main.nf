process GATK4_SVCLUSTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(vcfs), path(indices)
    path ploidy_table
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz"), emit: clustered_vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: clustered_vcf_index
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input = vcfs.collect { vcf -> "--variant ${vcf}" }.join(" ")

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK SVCluster] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SVCluster \\
        --output ${prefix}.vcf.gz \\
        --ploidy-table ${ploidy_table} \\
        ${input} \\
        --reference ${fasta} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
