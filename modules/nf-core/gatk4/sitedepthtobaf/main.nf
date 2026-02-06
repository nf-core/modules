process GATK4_SITEDEPTHTOBAF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(site_depths), path(site_depths_indices)
    tuple path(vcf), path(tbi)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.baf.txt.gz"), emit: baf
    tuple val(meta), path("*.baf.txt.gz.tbi"), emit: baf_tbi
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def site_depth_input = site_depths.collect { site_depth -> "--site-depth ${site_depth}" }.join(" ")
    def reference = fasta ? "--reference ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK SiteDepthtoBAF] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SiteDepthtoBAF \\
        --baf-evidence-output ${prefix}.baf.txt.gz \\
        --baf-sites-vcf ${vcf} \\
        ${site_depth_input} \\
        ${reference} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip -c > ${prefix}.baf.txt.gz
    touch ${prefix}.baf.txt.gz.tbi
    """
}
