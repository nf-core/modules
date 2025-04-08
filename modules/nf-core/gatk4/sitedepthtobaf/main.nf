process GATK4_SITEDEPTHTOBAF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(site_depths), path(site_depths_indices)
    tuple path(vcf), path(tbi)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.baf.txt.gz")       , emit: baf
    tuple val(meta), path("*.baf.txt.gz.tbi")   , emit: baf_tbi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def site_depth_input = site_depths.collect({"--site-depth ${it}"}).join(" ")
    def reference = fasta ? "--reference ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK SiteDepthtoBAF] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip -c > ${prefix}.baf.txt.gz
    touch ${prefix}.baf.txt.gz.tbi

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
