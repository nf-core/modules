process GATK4_SITEDEPTHTOBAF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(site_depths)
    path vcf
    path fasta
    path dict

    output:
    tuple val(meta), path("*.baf.txt.gz"), emit: baf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def site_depth_input = site_depths.collect({"--site_depth ${it}"}).join(" ")
    def reference = fasta ? "--reference ${fasta}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK SiteDepthtoBAF] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk --java-options "-Xmx${avail_mem}g" SiteDepthtoBAF \\
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
