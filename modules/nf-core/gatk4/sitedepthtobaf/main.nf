process GATK4_SITEDEPTHTOBAF {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

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
}
