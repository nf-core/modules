process GATK4_COLLECTSVEVIDENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(input), path(input_index), path(site_depth_vcf), path(site_depth_vcf_tbi)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.sr.txt.gz"), emit: split_read_evidence
    tuple val(meta), path("*.sr.txt.gz.tbi"), emit: split_read_evidence_index
    tuple val(meta), path("*.pe.txt.gz"), emit: paired_end_evidence
    tuple val(meta), path("*.pe.txt.gz.tbi"), emit: paired_end_evidence_index
    tuple val(meta), path("*.sd.txt.gz"), emit: site_depths, optional: true
    tuple val(meta), path("*.sd.txt.gz.tbi"), emit: site_depths_index, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def sd_vcf = site_depth_vcf ? "--sd-file ${prefix}.sd.txt.gz --site-depth-locs-vcf ${site_depth_vcf}" : ""
    def reference = fasta ? "--reference ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK COLLECTSVEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CollectSVEvidence \\
        ${args} \\
        --input ${input} \\
        --sr-file ${prefix}.sr.txt.gz \\
        --pe-file ${prefix}.pe.txt.gz \\
        ${sd_vcf} \\
        ${reference} \\
        --tmp-dir . \\
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sd_vcf = site_depth_vcf ? "echo '' | gzip > ${prefix}.sd.txt.gz" : ""
    def sd_vcf_tbi = site_depth_vcf_tbi ? "touch ${prefix}.sd.txt.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.sr.txt.gz
    touch ${prefix}.sr.txt.gz.tbi
    echo "" | gzip > ${prefix}.pe.txt.gz
    touch ${prefix}.pe.txt.gz.tbi
    ${sd_vcf}
    ${sd_vcf_tbi}
    """
}
